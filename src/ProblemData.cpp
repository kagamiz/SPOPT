#include "ProblemData.hpp"

namespace SPOPT {
    const double EPS = 1e-15;

    ProblemData::ProblemData(std::string fileName)
    {
        LoadConfig(fileName);
    }

    void ProblemData::LoadConfig(std::string fileName)
    {
        YAML::Node problemDataConfig = YAML::LoadFile(fileName);

        enableScaling            = problemDataConfig["enableScaling"].as<bool>(true);
        enableGradientConstraint = problemDataConfig["enableGradientConstraint"].as<bool>(false);

        objectiveFunction.LoadFromFile(problemDataConfig["objectiveFunctionFile"].as<std::string>("examples/affine.txt"));

        if (problemDataConfig["equalityConstraintFiles"]) {
            YAML::Node equalityConstraintFiles = problemDataConfig["equalityConstraintFiles"];
            for (int i = 0; i < equalityConstraintFiles.size(); i++) {
                originalEqualityConstraints.emplace_back(Polynomial(equalityConstraintFiles[i].as<std::string>()));
            }
        }

        if (problemDataConfig["inequalityConstraintFiles"]) {
            YAML::Node inequalityConstraintFiles = problemDataConfig["inequalityConstraintFiles"];
            for (int i = 0; i < inequalityConstraintFiles.size(); i++) {
                originalInequalityConstraints.emplace_back(Polynomial(inequalityConstraintFiles[i].as<std::string>()));
            }
        }

        if (problemDataConfig["IndexSetsFile"]) {
            std::ifstream ifs;
            ifs.open(problemDataConfig.as<std::string>("examples/index_sets.txt"));

            int groupNum;
            ifs >> groupNum;

            for (int i = 0; i < groupNum; i++) {
                int sz;
                ifs >> sz;
                std::vector<int> tmp(sz);
                for (int j = 0; j < sz; j++) {
                    ifs >> tmp[j]; tmp[j]--;
                }
                std::sort(tmp.begin(), tmp.end());
                originalIndexSets.emplace_back(tmp);
            }

            originalJunctionTree.resize(groupNum);
            for (int i = 0; i < groupNum - 1; i++) {
                int u, v;
                ifs >> u >> v;
                originalJunctionTree[u - 1].emplace_back(v - 1);
                originalJunctionTree[v - 1].emplace_back(u - 1);
            }
        }
    }

    void ProblemData::ConstructSDP(int degree)
    {
        hierarchyDegree = degree;

        _ConstructNewConstraints();
        _ConstructTermMap();
        _ConstructVectorB();
        _ConstructMatrixA();
        _ConstructVectorC();
        _ConstructScalingData();
    }

    void ProblemData::_ConstructNewConstraints()
    {
        std::vector<Polynomial> constraints;
        int eqNum = originalEqualityConstraints.size();
        int ineqNum = originalInequalityConstraints.size();

        constraints.insert(constraints.end(), originalEqualityConstraints.begin(), originalEqualityConstraints.end());
        constraints.insert(constraints.end(), originalInequalityConstraints.begin(), originalInequalityConstraints.end());

        int maxDegree = objectiveFunction.degree;

        for (auto poly : constraints) {
            maxDegree = std::max(maxDegree, poly.degree);
        }

        // construct a map for each term to determine which
        // index sets they belong
        std::map<std::vector<int>, int> indexSetToID;

        indexSetToID[{}] = 0;
        int ID = 0;
        for (auto originalIndexSet : originalIndexSets) {
            int indSetSize = originalIndexSet.size();
            for (int sz = 1; sz <= maxDegree; sz++) {
                // enumeration of all the subset of size
                // which is equal to `sz`
                for (int x = (1 << sz) - 1; x < (1 << indSetSize); ) {
                    std::vector<int> v;
                    for (int i = 0; i < indSetSize; i++) {
                        if ((x >> i) & 1) {
                            v.push_back(originalIndexSet[i]);
                        }
                    }
                    if (indexSetToID.find(v) == indexSetToID.end()) {
                        indexSetToID[v] = ID;
                    }
                    int t = x | (x - 1);
                    x = (t + 1) | (((~ t & - ~ t) - 1) >> (__builtin_ctz(x) + 1));
                }
            }
            ID++;
        }

        std::vector<std::vector<Monomial>> objectiveMonomials;
        std::vector<std::vector<std::vector<Monomial>>> constraintMonomials;

        std::vector<std::vector<int>> objectiveIDs, constraintIDs;

        int newVariableID = objectiveFunction.maxIndex + 1;

        objectiveMonomials.resize(originalIndexSets.size());
        for (auto monomial : objectiveFunction.monomials) {
            Term t = monomial.first;
            t.erase(unique(t.begin(), t.end()), t.end());
            objectiveMonomials[indexSetToID[t]].emplace_back(monomial.first, monomial.second, /* sorted = */true);
        }

        objectiveIDs.resize(originalIndexSets.size());
        for (int i = 0; i < originalIndexSets.size(); i++) {
            objectiveIDs[i].resize(originalIndexSets[i].size());
            for (int j = 0; j < originalIndexSets[i].size(); j++) {
                objectiveIDs[i][j] = newVariableID++;
            }
        }

        constraintMonomials.resize(constraints.size());
        constraintIDs.resize(constraints.size());

        for (int i = 0; i < constraints.size(); i++) {
            constraintMonomials[i].resize(originalIndexSets.size());
            constraintIDs[i].resize(originalIndexSets.size());
            for (auto monomial : constraints[i].monomials) {
                Term t = monomial.first;
                t.erase(unique(t.begin(), t.end()), t.end());
                constraintMonomials[i][indexSetToID[t]].emplace_back(monomial.first, monomial.second, /* sorted = */true);
            }

            for (int j = 0; j < originalIndexSets.size(); j++) {
                constraintIDs[i][j] = newVariableID++;
            }
        }

        std::vector<bool> visited(objectiveFunction.maxIndex + 1);
        int variableOrder = 1;
        _TraverseTree(0, -1, newVariableID, objectiveMonomials, objectiveIDs, constraintMonomials, constraintIDs, visited, variableOrder);
    }

    void ProblemData::_TraverseTree(int v, int p, int &newVariableIndex,
                                    std::vector<std::vector<Monomial>> &objectiveMonomials,
                                    std::vector<std::vector<int>> &objectiveIDs,
                                    std::vector<std::vector<std::vector<Monomial>>> &constraintMonomials,
                                    std::vector<std::vector<int>> &constraintIDs,
                                    std::vector<bool> &visited,
                                    int &variableOrder)
    {
        std::vector<int> newVariables, originalVariables;

        for (auto tmpID : objectiveIDs[v]) {
            if (variableOrderMap.find(tmpID) == variableOrderMap.end()) {
                variableOrderMap[tmpID] = -1;
                newVariables.emplace_back(tmpID);
            }
        }

        for (int i = 0; i < constraintIDs.size(); i++) {
            if (variableOrderMap.find(constraintIDs[i][v]) == variableOrderMap.end()) {
                variableOrderMap[constraintIDs[i][v]] = -1;
                newVariables.emplace_back(constraintIDs[i][v]);
            }
        }

        std::vector<int> newIndexSet = originalIndexSets[v];
        newIndexSet.insert(newIndexSet.end(), newVariables.begin(), newVariables.end());
        int newIndexSetID = convertedIndexSets.size();
        std::sort(newIndexSet.begin(), newIndexSet.end());
        convertedIndexSets.emplace_back(newIndexSet);

        if (p == -1) {
            for (int i = 0; i < constraintIDs.size(); i++) {
                Term t = {constraintIDs[i][v]};
                if (i < originalEqualityConstraints.size()) {
                    convertedEqualityConstraints.emplace_back(Polynomial(Monomial(t, /* sorted = */true)));
                    groupIDOfConvertedEqualityConstraints.emplace_back(newIndexSetID);
                }
                else {
                    convertedInequalityConstraints.emplace_back(Polynomial(Monomial(t, /* sorted = */true)));
                    groupIDOfConvertedInequalityConstraints.emplace_back(newIndexSetID);
                }
            }
        }

        for (int i = 0; i < originalIndexSets[v].size(); i++) {
            if (!visited[originalIndexSets[v][i]]) {
                visited[originalIndexSets[v][i]] = true;
                Term t = {objectiveIDs[v][i]};
                convertedEqualityConstraints.emplace_back(Polynomial(Monomial(t, /* sorted = */true)));
                groupIDOfConvertedEqualityConstraints.emplace_back(newIndexSetID);

                originalVariables.emplace_back(originalIndexSets[v][i]);
            }
        }

        int childNum = (int)originalJunctionTree[v].size() - (p == -1 ? 0 : 1);

        if (childNum == 0) {
            for (int i = 0; i < constraintMonomials.size(); i++) {
                Term t = {constraintIDs[i][v]};
                Polynomial poly(Monomial(t, -1, /* sorted = */true));
                for (int j = 0; j < constraintMonomials[i][v].size(); j++) {
                    poly += constraintMonomials[i][v][j];
                }
                convertedEqualityConstraints.emplace_back(poly);
                groupIDOfConvertedEqualityConstraints.emplace_back(newIndexSetID);
            }

            for (int i = 0; i < originalIndexSets[v].size(); i++) {
                Term t = {objectiveIDs[v][i]};
                Polynomial poly(Monomial(t, -1, /* sorted = */true));
                for (int j = 0; j < objectiveMonomials[v].size(); j++) {
                    Monomial dif = objectiveMonomials[v][j].DifferentiateBy(originalIndexSets[v][i]);
                    if (dif.term.size() == 0 && std::abs(dif.coefficient) <= EPS) continue;
                    poly += dif;
                }
                convertedEqualityConstraints.emplace_back(poly);
                groupIDOfConvertedEqualityConstraints.emplace_back(newIndexSetID);
            }
        }
        else {
            std::vector<int> prvConstraintVariables;
            for (int i = 0; i < constraintIDs.size(); i++) {
                prvConstraintVariables.emplace_back(constraintIDs[i][v]);
            }

            std::vector<int> prvObjectiveVariables = objectiveIDs[v];

            for (int i = 0; i < originalJunctionTree[v].size(); i++) {
                int nx = originalJunctionTree[v][i];
                if (nx == p) continue;

                std::vector<int> newIndexSet = originalIndexSets[v];

                newIndexSet.insert(newIndexSet.end(), prvConstraintVariables.begin(), prvConstraintVariables.end());
                newIndexSet.insert(newIndexSet.end(), prvObjectiveVariables.begin(), prvObjectiveVariables.end());

                for (int j = 0; j < constraintIDs.size(); j++) {
                    newIndexSet.emplace_back(constraintIDs[j][nx]);
                }

                int newIndexSetID = convertedIndexSets.size();

                std::vector<int> nxtConstraintVariables;

                for (int j = 0; j < constraintIDs.size(); j++) {
                    if (variableOrderMap.find(constraintIDs[j][nx]) == variableOrderMap.end()) {
                        variableOrderMap[constraintIDs[j][nx]] = -1;
                        newVariables.emplace_back(constraintIDs[j][nx]);
                    }

                    Polynomial poly;
                    if (childNum > 1) {
                        poly += Monomial(Term({prvConstraintVariables[j]}), -1, /* sorted = */true);
                        poly += Monomial(Term({constraintIDs[j][nx]}), 1, /* sorted = */true);
                        poly += Monomial(Term({newVariableIndex}), 1, /* sorted = */true);

                        newIndexSet.emplace_back(newVariableIndex);
                        nxtConstraintVariables.emplace_back(newVariableIndex);

                        variableOrderMap[newVariableIndex] = -1;
                        newVariables.emplace_back(newVariableIndex);
                        newVariableIndex++;
                    }
                    else {
                        poly += Monomial(Term({prvConstraintVariables[j]}), -1, /* sorted = */true);
                        poly += Monomial(Term({constraintIDs[j][nx]}), 1, /* sorted = */true);
                        for (int k = 0; k < constraintMonomials[j][v].size(); k++) {
                            poly += constraintMonomials[j][v][k];
                        }
                    }
                    convertedEqualityConstraints.emplace_back(poly);
                    groupIDOfConvertedEqualityConstraints.emplace_back(newIndexSetID);
                }

                std::vector<int> nxtObjectiveVariables;

                for (int j = 0; j < originalIndexSets[v].size(); j++) {
                    auto it = std::find(originalIndexSets[nx].begin(), originalIndexSets[nx].end(), originalIndexSets[v][j]);

                    if (it != originalIndexSets[nx].end()) {
                        int pos = it - originalIndexSets[nx].begin();

                        if (variableOrderMap.find(objectiveIDs[nx][pos]) == variableOrderMap.end()) {
                            variableOrderMap[objectiveIDs[nx][pos]] = -1;
                            newVariables.emplace_back(objectiveIDs[nx][pos]);
                        }

                        newIndexSet.emplace_back(objectiveIDs[nx][pos]);
                    }

                    Polynomial poly;

                    if (childNum > 1) {
                        if (it != originalIndexSets[nx].end()) {
                            int pos = it - originalIndexSets[nx].begin();

                            poly += Monomial(Term({prvObjectiveVariables[j]}), -1, /* sorted = */true);
                            poly += Monomial(Term({objectiveIDs[nx][pos]}), 1, /* sorted = */true);
                            poly += Monomial(Term({newVariableIndex}), 1, /* sorted = */true);

                            newIndexSet.emplace_back(newVariableIndex);
                            nxtObjectiveVariables.emplace_back(newVariableIndex);

                            variableOrderMap[newVariableIndex] = -1;
                            newVariables.emplace_back(newVariableIndex);
                            newVariableIndex++;
                        }
                        else {
                            nxtObjectiveVariables.emplace_back(prvObjectiveVariables[j]);
                        }
                    }
                    else {
                        poly += Monomial(Term({prvObjectiveVariables[j]}), -1, /* sorted = */true);
                        if (it != originalIndexSets[nx].end()) {
                            int pos = it - originalIndexSets[nx].begin();
                            poly += Monomial(Term({objectiveIDs[nx][pos]}), 1, /* sorted = */true);
                        }
                        for (int k = 0; k < objectiveMonomials[v].size(); k++) {
                            Monomial dif = objectiveMonomials[v][k].DifferentiateBy(originalIndexSets[v][j]);
                            if (dif.term.size() == 0 && std::abs(dif.coefficient) <= EPS) continue;
                            poly += dif;
                        }
                    }
                    if (poly.monomials.size() > 0) {
                        convertedEqualityConstraints.emplace_back(poly);
                        groupIDOfConvertedEqualityConstraints.emplace_back(newIndexSetID);
                    }
                }

                prvConstraintVariables = nxtConstraintVariables;
                prvObjectiveVariables  = nxtObjectiveVariables;

                childNum--;
                std::sort(newIndexSet.begin(), newIndexSet.end());
                convertedIndexSets.emplace_back(newIndexSet);
            }

            for (int i = 0; i < originalJunctionTree[v].size(); i++) {
                int nx = originalJunctionTree[v][i];
                if (nx == p) continue;
                _TraverseTree(nx, v, newVariableIndex, objectiveMonomials, objectiveIDs, constraintMonomials, constraintIDs, visited, variableOrder);
            }
        }

        std::reverse(newVariables.begin(), newVariables.end());

        for (int i = 0; i < newVariables.size(); i++) {
            variableOrderMap[newVariables[i]] = variableOrder++;
        }

        for (int i = 0; i < originalVariables.size(); i++) {
            variableOrderMap[originalVariables[i]] = variableOrder++;
        }
    }

    void ProblemData::_GenerateTerms(std::vector<int> &inds, std::vector<Term> &terms, Term &tmp, int pos, int left)
    {
        if (pos == inds.size()) {
            terms.push_back(tmp);
            return;
        }

        for (int i = 0; i <= left; i++) {
            for (int j = 0; j < i; j++) {
                tmp.push_back(inds[pos]);
            }
            _GenerateTerms(inds, terms, tmp, pos + 1, left - i);
            for (int j = 0; j < i; j++) {
                tmp.pop_back();
            }
        }
    }


    void ProblemData::_ConstructTermMap()
    {
        convertedTermSets.resize(convertedIndexSets.size());
        for (int i = 0; i < convertedIndexSets.size(); i++) {
            Term tmp;
            _GenerateTerms(convertedIndexSets[i], convertedTermSets[i], tmp, 0, hierarchyDegree);
        }

        std::vector<std::pair<int, Term>> terms;

        int counter = 0;
        for (int i = 0; i < convertedTermSets.size(); i++) {
            int sz = convertedTermSets[i].size();
            for (int j = 0; j < sz; j++) {
                for (int k = j; k < sz; k++) {
                    Term t1 = convertedTermSets[i][j];
                    Term t2 = convertedTermSets[i][k];
                    Term merged;
                    int it1 = 0, it2 = 0, len = t1.size() + t2.size();
                    int rank = 0;
                    while (it1 + it2 < len) {
                        if (it2 == t2.size() || (it1 < t1.size() && t1[it1] <= t2[it2])) {
                            rank = std::max(rank, variableOrderMap[t1[it1]]);
                            merged.push_back(t1[it1]); it1++;
                        }
                        else {
                            rank = std::max(rank, variableOrderMap[t2[it2]]);
                            merged.push_back(t2[it2]); it2++;
                        }
                    }
                    terms.emplace_back(rank, merged);
                }
            }
        }

        std::sort(terms.begin(), terms.end());
        terms.erase(std::unique(terms.begin(), terms.end()), terms.end());

        for (int i = 0; i < terms.size(); i++) {
            termToInteger[terms[i].second] = i;
        }
    }

    void ProblemData::_ConstructVectorB()
    {
        b.setZero(termToInteger.size());
        for (auto monomial : objectiveFunction.monomials) {
            b(termToInteger[monomial.first]) = monomial.second;
        }
    }

    void ProblemData::_ConstructMatrixA()
    {
        int rowNumOfA = termToInteger.size();
        int colNumOfA = 0;

        int maxSize = 0;
        for (auto convertedIndexSet : convertedIndexSets) {
            maxSize = std::max(maxSize, (int)convertedIndexSets.size());
        }

        std::vector<std::vector<long long>> C(maxSize + hierarchyDegree + 1, std::vector<long long>(maxSize + hierarchyDegree + 1));

        C[0][0] = 1;

        for (int i = 1; i <= maxSize + hierarchyDegree; i++) {
            C[i][0] = C[i][i] = 1;
            for (int j = 1; j < i; j++) {
                C[i][j] = C[i - 1][j] + C[i - 1][j - 1];
            }
        }

        std::vector<Polynomial> constraints;
        std::vector<int> groupIDs;
        
        constraints.insert(constraints.end(), convertedInequalityConstraints.begin(), convertedInequalityConstraints.end());
        constraints.insert(constraints.end(), convertedEqualityConstraints.begin(), convertedEqualityConstraints.end());

        groupIDs.insert(groupIDs.end(), groupIDOfConvertedInequalityConstraints.begin(), groupIDOfConvertedInequalityConstraints.end());
        groupIDs.insert(groupIDs.end(), groupIDOfConvertedEqualityConstraints.begin(), groupIDOfConvertedEqualityConstraints.end());
       
        colNumOfA += 1; // colNum of A_1
        // add colNum of A_2
        for (auto &convertedTermSet : convertedTermSets) {
            int sz = convertedTermSet.size();
            colNumOfA += sz * sz;
            psdMatrixSizes.emplace_back(sz);
        }
        // add colNum of A_3
        for (int i = 0; i < constraints.size(); i++) {
            int sz = convertedIndexSets[groupIDs[i]].size();
            int polyDeg = constraints[i].degree;
            int freeDeg = hierarchyDegree - (polyDeg + 1) / 2;
            colNumOfA += C[sz + freeDeg][freeDeg] * C[sz + freeDeg][freeDeg];

            if (i < convertedEqualityConstraints.size()) {
                psdMatrixSizes.emplace_back(C[sz + freeDeg][freeDeg]);
            }
            else {
                symmetricMatrixSizes.emplace_back(C[sz + freeDeg][freeDeg]);
            }
        }

        A.resize(rowNumOfA, colNumOfA);

        std::vector<Eigen::Triplet<double>> tripletList;
        tripletList.emplace_back(0, 0, 1);

        int leftmostPosition = 1;

        for (int i = 0; i < convertedTermSets.size(); i++) {
            int sz = convertedTermSets[i].size();
            for (int j = 0; j < sz; j++) {
                for (int k = 0; k < sz; k++) {
                    Term t1 = convertedTermSets[i][j];
                    Term t2 = convertedTermSets[i][k];
                    Term merged;
                    int it1 = 0, it2 = 0, len = t1.size() + t2.size();
                    while (it1 + it2 < len) {
                        if (it2 == t2.size() || (it1 < t1.size() && t1[it1] <= t2[it2])) {
                            merged.push_back(t1[it1]); it1++;
                        }
                        else {
                            merged.push_back(t2[it2]); it2++;
                        }
                    }
                    int id = termToInteger[merged];
                    tripletList.emplace_back(id, leftmostPosition, 1);
                    leftmostPosition++;
                }
            }
        }

        for (int pt = 0; pt < constraints.size(); pt++) {
            int i = groupIDs[pt];
            int sz = convertedTermSets[i].size();

            for (int j = 0; j < sz; j++) {
                for (int k = j; k < sz; k++) {
                    Term t1 = convertedTermSets[i][j];
                    Term t2 = convertedTermSets[i][k];
                    if (t1.size() + t2.size() + constraints[pt].degree > 2 * hierarchyDegree) continue;

                    Term _merged;
                    int it1 = 0, it2 = 0, len = t1.size() + t2.size();
                    while (it1 + it2 < len) {
                        if (it2 == t2.size() || (it1 < t1.size() && t1[it1] <= t2[it2])) {
                            _merged.push_back(t1[it1]); it1++;
                        }
                        else {
                            _merged.push_back(t2[it2]); it2++;
                        }
                    }
                    for (auto &monomial : constraints[pt].monomials) {
                        Term t3 = monomial.first;
                        double coef = monomial.second;
                        Term merged;
                        int it = 0, it3 = 0, len = _merged.size() + t3.size();

                        while (it + it3 < len) {
                            if (it3 == t3.size() || (it < _merged.size() && _merged[it] <= t3[it3])) {
                                merged.push_back(_merged[it]); it++;
                            }
                            else {
                                merged.push_back(t3[it3]); it3++;
                            }
                        }

                        int id = termToInteger[merged];
                        tripletList.emplace_back(id, leftmostPosition, coef);
                    }
                    leftmostPosition++;
                }
            }
        }
        A.setFromTriplets(tripletList.begin(), tripletList.end());
        //At = A.transpose();
        //AAt.compute(A * At);
    }

    void ProblemData::_ConstructVectorC()
    {
        c.resize(A.cols());
        c.coeffRef(0) = 1;
    }

    void ProblemData::_ConstructScalingData()
    {
        D.resize(A.rows());
        E.resize(A.cols());

        std::fill(D.begin(), D.end(), 1);
        std::fill(E.begin(), E.end(), 1);
        
        std::vector<double> Dt(A.rows()), Et(A.cols());

        const int iterationNum = 10;
        const double scalingMin = 1e-4;
        const double scalingMax = 1e+4;
        
        for (int repeat = 0; repeat < iterationNum; repeat++) {
            std::fill(Dt.begin(), Dt.end(), 0);
            std::fill(Et.begin(), Et.end(), 0);

            // calculate row norms and col norms
            for (int i = 0; i < A.outerSize(); i++) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
                    Dt[it.row()] = std::max(Dt[it.row()], std::abs(it.value()));
                    Et[it.col()] = std::max(Et[it.col()], std::abs(it.value()));
                }
            }

            for (int i = 0; i < A.rows(); i++) {
                Dt[i] = std::clamp(sqrt(Dt[i]), scalingMin / 10, scalingMax);
                if (Dt[i] < scalingMin) Dt[i] = 1.0;
            }

            for (int i = 0; i < A.cols(); i++) {
                Et[i] = std::clamp(sqrt(Et[i]), scalingMin / 10, scalingMax);
                if (Et[i] < scalingMin) Et[i] = 1.0;
            }

            // mean of E across each cone
            int leftmostPosition = 1;
            for (auto matSize : psdMatrixSizes) {
                double sum = 0;
                for (int i = 0; i < matSize * matSize; i++) {
                    sum += Et[leftmostPosition + i];
                }
                sum /= (matSize * matSize);

                for (int i = 0; i < matSize * matSize; i++) {
                    Et[leftmostPosition + i] = sum;
                }

                leftmostPosition += matSize * matSize;
            }

            // scale the rows and cols with Dt and Et
            for (int i = 0; i < A.outerSize(); i++) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
                    it.valueRef() /= (Dt[it.row()] * Et[it.col()]);
                }
            }

            // Accumulate scaling
            for (int i = 0; i < A.rows(); i++) {
                D[i] *= Dt[i];
            }
            for (int i = 0; i < A.cols(); i++) {
                E[i] *= Et[i];
            }
        }

        for (int i = 0; i < A.rows(); i++) {
            D[i] = 1 / D[i];
        }

        for (int i = 0; i < A.cols(); i++) {
            E[i] = 1 / E[i];
        }

        /*
        At = A.transpose();
        Eigen::SparseMatrix<double> AAtFull = A * At;
        AAt.compute(AAtFull);

        for (int i = 0; i < rowNumOfA; i++) {
            AAtFull.coeffRef(i, i) += 1;
        }
        IAAt.compute(AAtFull);
        */
        
        // calculate mean of row / col norms of A
        std::vector<double> rowNorms(A.rows(), 0), colNorms(A.cols(), 0);
        for (int i = 0; i < A.outerSize(); i++) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
                rowNorms[it.row()] += it.value() * it.value();
                colNorms[it.col()] += it.value() * it.value();
            }
        }

        double rowNormMean = 0, colNormMean = 0;

        for (int i = 0; i < A.rows(); i++) {
            rowNormMean += sqrt(rowNorms[i]) / A.rows();
        }

        for (int i = 0; i < A.cols(); i++) {
            colNormMean += sqrt(colNorms[i]) / A.cols();
        }

        // scale f
        for (int i = 0; i < A.rows(); i++) {
            b(i) *= D[i];
        }
        dualScaler = colNormMean / std::max(b.norm(), 1e-6);
        for (int i = 0; i < A.rows(); i++) {
            b(i) *= dualScaler;
        }

        // scale c
        c *= E[0];
        primalScaler = rowNormMean / std::max(c.norm(), 1e-6);
        c *= primalScaler;
    }
}