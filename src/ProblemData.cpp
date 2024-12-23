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

        LoadConfig(problemDataConfig);
    }

    void ProblemData::LoadConfig(YAML::Node &problemDataConfig)
    {
        enableScaling                   = problemDataConfig["enableScaling"].as<bool>(true);
        verbose                         = problemDataConfig["verbose"].as<bool>(false);
        scalingFactor                   = problemDataConfig["scalingFactor"].as<double>(5.0);
        gradientConstraintType          = problemDataConfig["gradientConstraintType"].as<int>(0);
        enableLowerBound                = problemDataConfig["enableLowerBound"].as<bool>(false);
        lowerBoundConstant              = problemDataConfig["lowerBoundConstant"].as<double>(0.0);
        enableUpperBound                = problemDataConfig["enableUpperBound"].as<bool>(false);
        upperBoundConstant              = problemDataConfig["upperBoundConstant"].as<double>(0.0);
        enableOriginalVariableNormBound = problemDataConfig["enableOriginalVariableNormBound"].as<bool>(false);
        originalVariableNormBound       = problemDataConfig["originalVariableNormBound"].as<double>(0.0);
        perturbObjectiveFunction        = problemDataConfig["perturbObjectiveFunction"].as<bool>(false);
        addVariableNonnegativity        = problemDataConfig["addVariableNonnegativity"].as<bool>(false);
        addFirstOrderFullMomentMatrix   = problemDataConfig["addFirstOrderFullMomentMatrix"].as<bool>(false);
        RNGSeed                         = problemDataConfig["RNGSeed"].as<unsigned int>(0);

        srand(RNGSeed);

        objectiveFunction.LoadFromFile(problemDataConfig["objectiveFunctionFile"].as<std::string>("examples/affine.txt"));

        if (problemDataConfig["maximize"].as<bool>(false) == true) {
            objectiveFunction = -objectiveFunction;
            std::swap(enableLowerBound, enableUpperBound);
            std::swap(lowerBoundConstant, upperBoundConstant);
            lowerBoundConstant *= -1;
            upperBoundConstant *= -1;
        }

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

        if (problemDataConfig["indexSetsFile"]) {
            std::ifstream ifs;
            ifs.open(problemDataConfig["indexSetsFile"].as<std::string>("examples/index_sets.txt"));

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
        else {
            originalJunctionTree.resize(1);
            originalJunctionTree[0].clear();
            int n = objectiveFunction.maxIndex + 1;
            std::vector<int> allIndex(n);
            std::iota(allIndex.begin(), allIndex.end(), 0);
            originalIndexSets.emplace_back(allIndex);
        }

        hierarchyDegree = problemDataConfig["hierarchyDegree"].as<int>((objectiveFunction.degree + 1) / 2);
    }

    bool ProblemData::IsUnconstrained() const
    {
        return    convertedEqualityConstraints.size() == 0 && convertedInequalityConstraints.size() == 0
               && gradientConstraintType == 0;
    }

    void ProblemData::ConstructSparsification()
    {
        if (gradientConstraintType == 2) {
            _AddGradientConstraints();
        }

        if (gradientConstraintType == 3) {
            _AddMinors();
        }
        _ConstructNewConstraints();
    }

    void ProblemData::ConstructSDP()
    {
        if (perturbObjectiveFunction) {
            int n = objectiveFunction.maxIndex + 1;
            Eigen::VectorXd r = Eigen::VectorXd::Random(n);
            r /= r.norm();
            r /= 10.0;
            for (int i = 0; i < n; i++) {
                objectiveFunction += Monomial(Term({i}), r(i), true);
            }
        }

        ConstructSparsification();

	int n = objectiveFunction.maxIndex;
        std::set<int> visited;
        for (int i = 0; i < convertedIndexSets.size(); i++) {
            auto &convertedIndexSet = convertedIndexSets[i];
            for (auto &e : convertedIndexSet) {
                if (visited.find(e) == visited.end()) {
                    visited.insert(e);
                    if (addVariableNonnegativity && e <= n) {
                        convertedInequalityConstraints.emplace_back(Term({e}));
                        groupIDOfConvertedInequalityConstraints.emplace_back(i);
                    }

                    if (enableOriginalVariableNormBound) {
                        convertedInequalityConstraints.emplace_back(Polynomial(Monomial(originalVariableNormBound * originalVariableNormBound)) - Monomial(Term({e, e})));
                        groupIDOfConvertedInequalityConstraints.emplace_back(i);
                    }
                }
            }
        }

        _ConstructTermMap();
        _ConstructMatrixA();
        _ConstructVectorB();
        _ConstructVectorC();
        _ConstructScalingData();

        std::map<int, int> PSDSize;
        for (auto psdMat : psdMatrixSizes) {
            if (PSDSize.find(psdMat) != PSDSize.end()) PSDSize[psdMat]++;
            else PSDSize[psdMat] = 1;
        }

        if (verbose) {
            std::cout << "semidefinite variables :" << std::endl;
            for (auto it = PSDSize.begin(); it != PSDSize.end(); it++) {
                auto [p, q] = *it;
                std::cout << "(" << p << "," << q << ")" << std::endl;
            }
            std::cout << "Number of the LMI : " << A.rows() << std::endl;
        }
    }

    void ProblemData::_AddGradientConstraints()
    {
        int variableNum = objectiveFunction.maxIndex + 1;
        Polynomial LagrangianFunction = objectiveFunction;
        
        for (int i = 0; i < variableNum; i++) {
            LagrangianFunction += Monomial(Term({i, i, variableNum}), 1, /* sorted = */true);
        }
        LagrangianFunction += Monomial(Term({variableNum}), -1, /* sorted = */true);

        variableNum++;
        for (int i = 0; i < variableNum; i++) {
            double D = 1;
            Polynomial f = LagrangianFunction.DifferentiateBy(i);
            for (auto &monomial : f.monomials) {
                D = std::max(std::abs(D), monomial.second);
            }
            for (auto &monomial : f.monomials) {
                monomial.second /= D;
            }
            originalEqualityConstraints.emplace_back(f);
        }

        for (int i = 0; i < originalIndexSets.size(); i++) {
            originalIndexSets[i].emplace_back(variableNum - 1);
        }
	    //Polynomial TotalNorm = Polynomial(Monomial(1)) - Monomial(Term({variableNum - 1, variableNum - 1}), 1./10000, /* sorted = */true);
        //originalInequalityConstraints.emplace_back(TotalNorm);
	    objectiveFunction.maxIndex++;
	    //objectiveFunction = LagrangianFunction;
    }

    void ProblemData::_AddMinors()
    {
        int n = objectiveFunction.maxIndex + 1;

        originalJunctionTree.resize(1);
        originalJunctionTree[0].clear();
        originalIndexSets.clear();
        std::vector<int> allIndex(n);
        std::iota(allIndex.begin(), allIndex.end(), 0);
        originalIndexSets.emplace_back(allIndex);

        std::vector<Polynomial> pderivatives;

        for (int i = 0; i < n; i++) {
            pderivatives.emplace_back(objectiveFunction.DifferentiateBy(i));
        }

        std::vector<Polynomial> newConstraints(2 * n - 3);
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                newConstraints[i + j - 1] += pderivatives[i] * Monomial(Term({j})) - pderivatives[j] * Monomial(Term({i}));
            }
        }

        for (int i = 0; i < 2 * n - 3; i++) {
            newConstraints[i].Simplify();
            double D = 1;
            for (auto &monomial : newConstraints[i].monomials) D = std::max(std::abs(D), monomial.second);
            for (auto &monomial : newConstraints[i].monomials) monomial.second /= D;
	    }

        originalEqualityConstraints.insert(originalEqualityConstraints.end(), newConstraints.begin(), newConstraints.end());
    }

    void ProblemData::_ConstructNewConstraints()
    {
        std::vector<Polynomial>            constraints;
        std::vector<ConstraintType>        constraintTypes;
        std::vector<Polynomial>            unpartitionedConstraints;
        std::vector<ConstraintType>        unpartitionedConstraintTypes;

        if (enableLowerBound) {
            Polynomial lb = objectiveFunction - Monomial(lowerBoundConstant);
            double D = 1;
            for (auto &monomial : lb.monomials) D = std::max(std::abs(monomial.second), D);
            for (auto &monomial : lb.monomials) monomial.second /= D;
            originalInequalityConstraints.emplace_back(lb);
        }
	
	    if (enableUpperBound) {
            Polynomial ub = -(objectiveFunction - Monomial(upperBoundConstant));
            double D = 1;
            for (auto &monomial : ub.monomials) D = std::max(std::abs(monomial.second), D);
            for (auto &monomial : ub.monomials) monomial.second /= D;
            originalInequalityConstraints.emplace_back(ub);
        }
	
        std::vector<Polynomial> tmpPolys;
        tmpPolys.insert(tmpPolys.end(), originalEqualityConstraints.begin(), originalEqualityConstraints.end());
        tmpPolys.insert(tmpPolys.end(), originalInequalityConstraints.begin(), originalInequalityConstraints.end());

        for (int i = 0; i < tmpPolys.size(); i++) {
            auto &originalConstraint = tmpPolys[i];
            std::set<int> variableIndices;
            for (auto &monomial : originalConstraint.monomials) {
                for (auto ind : monomial.first) {
                    variableIndices.insert(ind);
                }
            }
            bool isSubsetOfOriginalSet = false;
            for (auto &originalIndexSet : originalIndexSets) {
                std::set<int> indSet(originalIndexSet.begin(), originalIndexSet.end());
                if (std::includes(indSet.begin(), indSet.end(), variableIndices.begin(), variableIndices.end())) {
                    isSubsetOfOriginalSet = true;
                    break;
                }
            }

            if (isSubsetOfOriginalSet) {
                unpartitionedConstraints.emplace_back(originalConstraint);
                unpartitionedConstraintTypes.emplace_back(i < originalEqualityConstraints.size() ? ConstraintType::EqualityConstraint : ConstraintType::InequalityConstraint);
            }
            else {
                constraints.emplace_back(originalConstraint);
                constraintTypes.emplace_back(i < originalEqualityConstraints.size() ? ConstraintType::EqualityConstraint : ConstraintType::InequalityConstraint);
            }
        }

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

        if (gradientConstraintType == 1) {
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

            for (int j = 1; j < originalIndexSets.size(); j++) {
                constraintIDs[i][j] = newVariableID++;
            }
        }

        std::vector<bool> visited(objectiveFunction.maxIndex + 1);
        int variableOrder = 0;
        _TraverseTree(0, -1, newVariableID, objectiveMonomials, objectiveIDs, constraintMonomials, constraintIDs, constraintTypes, visited, variableOrder);

        for (int i = 0; i < unpartitionedConstraints.size(); i++) {
            auto &unpartitionedConstraint = unpartitionedConstraints[i];
            std::set<int> variableIndices;
            for (auto &monomial : unpartitionedConstraint.monomials) {
                for (auto ind : monomial.first) {
                    variableIndices.insert(ind);
                }
            }
            
            bool foundSet = false;
            for (int j = 0; j < convertedIndexSets.size(); j++) {
                std::set<int> indSet(convertedIndexSets[j].begin(), convertedIndexSets[j].end());
                if (std::includes(indSet.begin(), indSet.end(), variableIndices.begin(), variableIndices.end())) {
                    foundSet = true;
                    if (unpartitionedConstraintTypes[i] == ConstraintType::EqualityConstraint) {
                        convertedEqualityConstraints.emplace_back(unpartitionedConstraint);
                        groupIDOfConvertedEqualityConstraints.emplace_back(j);
                    }
                    else {
                        convertedInequalityConstraints.emplace_back(unpartitionedConstraint);
                        groupIDOfConvertedInequalityConstraints.emplace_back(j);
                    }
		            break;
                }
            }
            if (!foundSet) {
                std::cout << "[ERROR] preprocessing of polynomial decomposition failed!!" << std::endl;
                exit(EXIT_FAILURE);
            }
        }

    }

    void ProblemData::_TraverseTree(int v, int p, int &newVariableIndex,
                                    std::vector<std::vector<Monomial>> &objectiveMonomials,
                                    std::vector<std::vector<int>> &objectiveIDs,
                                    std::vector<std::vector<std::vector<Monomial>>> &constraintMonomials,
                                    std::vector<std::vector<int>> &constraintIDs,
                                    std::vector<ConstraintType>  &constraintTypes,
                                    std::vector<bool> &visited,
                                    int &variableOrder)
    {
        std::vector<int> newVariables, originalVariables;

        if (gradientConstraintType == 1) {
            for (auto tmpID : objectiveIDs[v]) {
                if (variableOrderMap.find(tmpID) == variableOrderMap.end()) {
                    variableOrderMap[tmpID] = -1;
                    newVariables.emplace_back(tmpID);
                }
            }
        }

        if (v != 0) {
            for (int i = 0; i < constraintIDs.size(); i++) {
                if (variableOrderMap.find(constraintIDs[i][v]) == variableOrderMap.end()) {
                    variableOrderMap[constraintIDs[i][v]] = -1;
                    newVariables.emplace_back(constraintIDs[i][v]);
                }
            }
        }

        int newIndexSetID = convertedIndexSets.size();
        #if 0
        if (p == -1) {
            for (int i = 0; i < constraintIDs.size(); i++) {
                Term t = {constraintIDs[i][v]};
                if (constraintTypes[i] == ConstraintType::EqualityConstraint) {
                    convertedEqualityConstraints.emplace_back(Polynomial(Monomial(t, /* sorted = */true)));
                    groupIDOfConvertedEqualityConstraints.emplace_back(newIndexSetID);
                }
                else if (constraintTypes[i] == ConstraintType::InequalityConstraint) {
                    convertedInequalityConstraints.emplace_back(Polynomial(Monomial(t, /* sorted = */true)));
                    groupIDOfConvertedInequalityConstraints.emplace_back(newIndexSetID);
                }
                else {
                    convertedInequalityConstraints.emplace_back(Polynomial(Monomial(t, /* sorted = */true)) - Polynomial(Monomial(lowerBoundConstant)));
                    groupIDOfConvertedInequalityConstraints.emplace_back(newIndexSetID);

                    convertedInequalityConstraints.emplace_back(Polynomial(Monomial(upperBoundConstant)) - Polynomial(Monomial(t, /* sorted = */true)));
                    groupIDOfConvertedInequalityConstraints.emplace_back(newIndexSetID);
                }
            }
        }

        for (int i = 0; i < originalIndexSets[v].size(); i++) {
            if (!visited[originalIndexSets[v][i]]) {
                visited[originalIndexSets[v][i]] = true;

                if (gradientConstraintType == 1) {
                    Term t = {objectiveIDs[v][i]};
                    convertedEqualityConstraints.emplace_back(Polynomial(Monomial(t, /* sorted = */true)));
                    groupIDOfConvertedEqualityConstraints.emplace_back(newIndexSetID);
                }

                originalVariables.emplace_back(originalIndexSets[v][i]);
            }
        }
        #endif

        for (int i = 0; i < originalIndexSets[v].size(); i++) {
            if (!visited[originalIndexSets[v][i]]) {
                visited[originalIndexSets[v][i]] = true;
                originalVariables.emplace_back(originalIndexSets[v][i]);
            }
        }

        int childNum = (int)originalJunctionTree[v].size() - (p == -1 ? 0 : 1);

        if (childNum == 0) {
            std::vector<int> newIndexSet = originalIndexSets[v];
            if (v != 0) {
                for (int i = 0; i < constraintIDs.size(); i++) {
                    newIndexSet.emplace_back(constraintIDs[i][v]);
                }
            }
            if (gradientConstraintType == 1) {
                newIndexSet.insert(newIndexSet.end(), objectiveIDs[v].begin(), objectiveIDs[v].end());
            }
            std::sort(newIndexSet.begin(), newIndexSet.end());
            convertedIndexSets.emplace_back(newIndexSet);

            for (int i = 0; i < constraintMonomials.size(); i++) {
                Polynomial poly;
                if (v != 0) {
                    Term t = {constraintIDs[i][v]};
                    poly += Monomial(t, -1, /* sorted = */true);
                }
                for (int j = 0; j < constraintMonomials[i][v].size(); j++) {
                    poly += constraintMonomials[i][v][j];
                }

                if (v == 0 && constraintTypes[i] == ConstraintType::InequalityConstraint) {
                    convertedInequalityConstraints.emplace_back(poly);
                    groupIDOfConvertedInequalityConstraints.emplace_back(newIndexSetID);
                }
                else {
                    convertedEqualityConstraints.emplace_back(poly);
                    groupIDOfConvertedEqualityConstraints.emplace_back(newIndexSetID);
                }
            }

            if (gradientConstraintType == 1) {
                for (int i = 0; i < originalIndexSets[v].size(); i++) {
                    Term t = {objectiveIDs[v][i]};
                    Polynomial poly(Monomial(t, -1, /* sorted = */true));
                    for (int j = 0; j < objectiveMonomials[v].size(); j++) {
                        Monomial dif = objectiveMonomials[v][j].DifferentiateBy(originalIndexSets[v][i]);
                        poly += dif;
                    }
                    poly.Simplify();
                    convertedEqualityConstraints.emplace_back(poly);
                    groupIDOfConvertedEqualityConstraints.emplace_back(newIndexSetID);
                }
            }
        }
        else {
            std::vector<int> prvConstraintVariables;
            if (v != 0) {
                for (int i = 0; i < constraintIDs.size(); i++) {
                    prvConstraintVariables.emplace_back(constraintIDs[i][v]);
                }
            }

            std::vector<int> prvObjectiveVariables = (gradientConstraintType == 1 ? objectiveIDs[v] : std::vector<int>());

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
                        if (!(v == 0 && i == 0)) {
                            poly += Monomial(Term({prvConstraintVariables[j]}), -1, /* sorted = */true);
                        }
                        poly += Monomial(Term({constraintIDs[j][nx]}), 1, /* sorted = */true);
                        poly += Monomial(Term({newVariableIndex}), 1, /* sorted = */true);

                        newIndexSet.emplace_back(newVariableIndex);
                        nxtConstraintVariables.emplace_back(newVariableIndex);

                        variableOrderMap[newVariableIndex] = -1;
                        newVariables.emplace_back(newVariableIndex);
                        newVariableIndex++;
                    }
                    else {
                        if (!(v == 0 && i == 0)) {
                            poly += Monomial(Term({prvConstraintVariables[j]}), -1, /* sorted = */true);
                        }
                        poly += Monomial(Term({constraintIDs[j][nx]}), 1, /* sorted = */true);
                        for (int k = 0; k < constraintMonomials[j][v].size(); k++) {
                            poly += constraintMonomials[j][v][k];
                        }
                    }
                    if (v == 0 && constraintTypes[j] == ConstraintType::InequalityConstraint) {
                        convertedInequalityConstraints.emplace_back(poly);
                        groupIDOfConvertedInequalityConstraints.emplace_back(newIndexSetID);
                    }
                    else {
                        convertedEqualityConstraints.emplace_back(poly);
                        groupIDOfConvertedEqualityConstraints.emplace_back(newIndexSetID);
                    }
                }

                std::vector<int> nxtObjectiveVariables;

                if (gradientConstraintType == 1) {
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
                                poly += dif;
                            }
                            poly.Simplify();
                        }
                        if (poly.monomials.size() > 0) {
                            convertedEqualityConstraints.emplace_back(poly);
                            groupIDOfConvertedEqualityConstraints.emplace_back(newIndexSetID);
                        }
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
                _TraverseTree(nx, v, newVariableIndex, objectiveMonomials, objectiveIDs, constraintMonomials, constraintIDs, constraintTypes, visited, variableOrder);
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
	
	//int firstAuxVarID = objectiveFunction.maxIndex + 1;
        for (int i = 0; i < convertedTermSets.size(); i++) {
            int sz = convertedTermSets[i].size();
            for (int j = 0; j < sz; j++) {
                for (int k = j; k < sz; k++) {
                    Term t1 = convertedTermSets[i][j];
                    Term t2 = convertedTermSets[i][k];
		    //if (t1.size() + t2.size() == 2 * hierarchyDegree && std::min(t1[0], t2[0]) >= firstAuxVarID) continue;
                    Term merged;
                    int it1 = 0, it2 = 0, len = t1.size() + t2.size();
                    int rank = variableOrderMap.size();
                    while (it1 + it2 < len) {
                        if (it2 == t2.size() || (it1 < t1.size() && t1[it1] <= t2[it2])) {
                            rank = std::min(rank, variableOrderMap[t1[it1]]);
                            merged.push_back(t1[it1]); it1++;
                        }
                        else {
                            rank = std::min(rank, variableOrderMap[t2[it2]]);
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

        if (addFirstOrderFullMomentMatrix) {
            int n = objectiveFunction.maxIndex + 1;
            int sz = termToInteger.size();
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    if (termToInteger.find({i, j}) == termToInteger.end()) {
                        termToInteger[{i, j}] = sz++;
                    }
                }
            }
        }

    }

    void ProblemData::_ConstructVectorB()
    {
        b.setZero(A.rows());

        for (auto monomial : objectiveFunction.monomials) {
            b(termToInteger[monomial.first]) = monomial.second;
        }

        originalBNorm = b.norm();
    }

    void ProblemData::_ConstructMatrixA()
    {
        int rowNumOfA = termToInteger.size();
        int colNumOfA = 0;

        int maxSize = 0;
        for (auto convertedIndexSet : convertedIndexSets) {
            maxSize = std::max(maxSize, (int)convertedIndexSet.size());
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
        if (addFirstOrderFullMomentMatrix) {
            int n = objectiveFunction.maxIndex + 1;
            colNumOfA += (n + 1) * (n + 2) / 2;
            psdMatrixSizes.emplace_back(n + 1);
        }
	//int firstAuxVarID = objectiveFunction.maxIndex + 1;
        for (auto &convertedTermSet : convertedTermSets) {
            //int sz = std::count_if(convertedTermSet.begin(), convertedTermSet.end(), [&](Term x) { return !(x.size() == hierarchyDegree && x[0] >= firstAuxVarID ); } );
            int sz = convertedTermSet.size();
	    colNumOfA += sz * (sz + 1) / 2;
            psdMatrixSizes.emplace_back(sz);
        }
        // add colNum of A_3
        for (int i = 0; i < constraints.size(); i++) {
            int sz = convertedIndexSets[groupIDs[i]].size();
            int polyDeg = constraints[i].degree;
            int freeDeg = hierarchyDegree - (polyDeg + 1) / 2;
            colNumOfA += C[sz + freeDeg][freeDeg] * (C[sz + freeDeg][freeDeg] + 1) / 2;

            if (i < convertedInequalityConstraints.size()) {
                psdMatrixSizes.emplace_back(C[sz + freeDeg][freeDeg]);
            }
            else {
                symmetricMatrixSizes.emplace_back(C[sz + freeDeg][freeDeg]);
            }
        }

        A.resize(rowNumOfA, colNumOfA);

        std::vector<Eigen::Triplet<double>> tripletList;

        tripletList.emplace_back(termToInteger[Term({})], 0, 1);

        int leftmostPosition = 1;
        double sq2 = sqrt(2);
        if (addFirstOrderFullMomentMatrix) {
            int n = objectiveFunction.maxIndex + 1;
            for (int i = 0; i <= n; i++) {
                for (int j = i; j <= n; j++) {
                    Term t = (i == 0 ? (j == 0 ? Term({}) : Term({j - 1})) : Term({i - 1, j - 1}));
                    int id = termToInteger[t];
                    tripletList.emplace_back(id, leftmostPosition, 1 * (i != j ? sq2 : 1));
                    if (i != j) sq2Cols.emplace_back(leftmostPosition);
                    leftmostPosition++;
                }
            }
        }
	
        for (int i = 0; i < convertedTermSets.size(); i++) {
            int sz = convertedTermSets[i].size();
            for (int j = 0; j < sz; j++) {
                for (int k = j; k < sz; k++) {
                    Term t1 = convertedTermSets[i][j];
                    Term t2 = convertedTermSets[i][k];
                    Term merged;
		    //if (t1.size() == hierarchyDegree && t1[0] >= firstAuxVarID) continue;
		    //if (t2.size() == hierarchyDegree && t2[0] >= firstAuxVarID) continue;
                    int it1 = 0, it2 = 0, len = t1.size() + t2.size();
                    while (it1 + it2 < len) {
                        if (it2 == t2.size() || (it1 < t1.size() && t1[it1] <= t2[it2])) {
                            merged.push_back(t1[it1]); it1++;
                        }
                        else {
                            merged.push_back(t2[it2]); it2++;
                        }
                    }
		    if (termToInteger.find(merged) == termToInteger.end()) {
		        for (auto &e : merged) {
			    std::cout << e << ", ";
			}
			std::cout << std::endl;
		    }
                    int id = termToInteger[merged];
                    tripletList.emplace_back(id, leftmostPosition, 1 * (j != k ? sq2 : 1));
                    if (j != k) sq2Cols.emplace_back(leftmostPosition);
                    leftmostPosition++;
                }
            }
        }

        for (int pt = 0; pt < constraints.size(); pt++) {
            int i = groupIDs[pt];
            int sz = convertedTermSets[i].size();
            int polyDeg = constraints[pt].degree;

            for (int j = 0; j < sz; j++) {
                for (int k = j; k < sz; k++) {
                    Term t1 = convertedTermSets[i][j];
                    Term t2 = convertedTermSets[i][k];
                    if (t1.size() > hierarchyDegree - (polyDeg + 1) / 2 || t2.size() > hierarchyDegree - (polyDeg + 1) / 2) continue;

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
                        if (termToInteger.find(merged) == termToInteger.end()) {
                            for (auto &e : merged) {
                                std::cout << e << ", ";
                            }
                            std::cout << std::endl;
                        }
                        int id = termToInteger[merged];
                        tripletList.emplace_back(id, leftmostPosition, coef * (j != k ? sq2 : 1));
                    }
                    if (j != k) sq2Cols.emplace_back(leftmostPosition);
                    leftmostPosition++;
                }
            }
        }
        A.setFromTriplets(tripletList.begin(), tripletList.end());
    }

    void ProblemData::_ConstructVectorC()
    {
        c.resize(A.cols());
        c.coeffRef(0) = 1;

        originalCNorm = 1;
    }

    void ProblemData::_ConstructScalingData()
    {
        D.resize(A.rows());
        E.resize(A.cols());

        std::fill(D.begin(), D.end(), 1);
        std::fill(E.begin(), E.end(), 1);

        if (enableScaling) {
            
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
                        //Dt[it.row()] += it.value() * it.value();
                        //Et[it.col()] += it.value() * it.value();
                    }
                }

                for (int i = 0; i < A.cols(); i++) {
                    Et[i] = std::clamp(sqrt(Et[i]), scalingMin / 10, scalingMax);
                    if (Et[i] < scalingMin) Et[i] = 1.0;
                }

                for (int i = 0; i < A.rows(); i++) {
                    Dt[i] = std::clamp(sqrt(Dt[i]), scalingMin / 10, scalingMax);
                    if (Dt[i] < scalingMin) Dt[i] = 1.0;
                }
                /*
                for (int i = 0; i < A.cols(); i++) {
                    Et[i] = Et[i] * A.cols() / A.rows();
                }
                */

                // mean of E across each cone
                int leftmostPosition = 1;
                for (auto matSize : psdMatrixSizes) {
                    double sum = 0;
                    for (int i = 0; i < matSize * (matSize + 1) / 2; i++) {
                        sum += Et[leftmostPosition + i];
                    }
                    sum /= (matSize * matSize);

                    for (int i = 0; i < matSize * (matSize + 1) / 2; i++) {
                        Et[leftmostPosition + i] = sum;
                    }

                    leftmostPosition += matSize * (matSize + 1) / 2;
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
        }
        
        primalScaler = 1;
        dualScaler = 1;

        if (enableScaling) {
            // calculate mean of row / col norms of A
            std::vector<double> rowNorms(A.rows(), 0), colNorms(A.cols(), 0);
            std::vector<double> rowInfNorms(A.rows(), 0), colInfNorms(A.cols(), 0);
            for (int i = 0; i < A.outerSize(); i++) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
                    rowNorms[it.row()] += it.value() * it.value();
                    rowInfNorms[it.row()] = std::max(rowInfNorms[it.row()], abs(it.value()));
                    colNorms[it.col()] += it.value() * it.value();
                    colInfNorms[it.col()] = std::max(colInfNorms[it.col()], abs(it.value()));
                }
            }

            double rowNormMean = 0, colNormMean = 0;

            for (int i = 0; i < A.rows(); i++) {
                rowNormMean += sqrt(rowNorms[i]) / A.rows();
                //rowNormMean += rowInfNorms[i] / A.rows();
            }

            for (int i = 0; i < A.cols(); i++) {
                colNormMean += sqrt(colNorms[i]) / A.cols();
                //colNormMean += colInfNorms[i] / A.cols();
            }

            for (int i = 0; i < A.outerSize(); i++) {
                for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
                    it.valueRef() *= scalingFactor;
                }
            }

            // scale b
            for (int i = 0; i < A.rows(); i++) {
                b(i) *= D[i];
            }
            dualScaler = colNormMean / std::max(b.norm(), 1e-6);
            b *= dualScaler * scalingFactor;

            // scale c
            for (Eigen::SparseVector<double>::InnerIterator it(c); it; ++it) {
                it.valueRef() *= E[it.index()];
            }

            primalScaler = rowNormMean / std::max(c.norm(), 1e-6);
            c *= primalScaler * scalingFactor;
        }
    }

    void ProblemData::ShowPOP()
    {
        std::cout << "minimize " << objectiveFunction.ToString() << std::endl;
        for (auto &eqConstraint : convertedEqualityConstraints) {
            std::cout << eqConstraint.ToString() << " == 0" << std::endl;
        }

        for (auto &ineqConstraint : convertedInequalityConstraints) {
            std::cout << ineqConstraint.ToString() << " >= 0" << std::endl;
        }
    }

    void ProblemData::OutputJuliaFile(std::string fileName)
    {
        std::ostream &os = (fileName != "" ? *(new std::ofstream(fileName)) : std::cout);

        os << "using CPUTime;" << std::endl;
        os << "using TSSOS;" << std::endl;
        os << "using DynamicPolynomials;" << std::endl;
        os << "using SparseArrays;" << std::endl;
        os << "using MultivariatePolynomials;" << std::endl;

        int maxIndex = objectiveFunction.maxIndex;
        for (auto convertedEqualityConstraint : convertedEqualityConstraints) {
            maxIndex = std::max(maxIndex, (int)convertedEqualityConstraint.maxIndex);
        }
        for (auto convertedInequalityConstraint : convertedInequalityConstraints) {
            maxIndex = std::max(maxIndex, (int)convertedInequalityConstraint.maxIndex);
        }

        os << "@polyvar x[1:" << maxIndex + 1 << "];" << std::endl;
        //os << "f=" << objectiveFunction.ToString(/* oneIndexed = */true) << ";" << std::endl;
        os << "f=x[1];" << std::endl;
        for (auto &monomial : objectiveFunction.monomials) {
            os << "f+=" << Monomial(monomial.first, monomial.second, /* sorted = */true).ToString(/* oneIndexed = */true, /* withSign = */false) << ";" << std::endl;
        }
        os << "f-=x[1];" << std::endl;
        os << "pop=[f];" << std::endl;
        for (auto convertedInequalityConstraint : convertedInequalityConstraints) {
            os << "push!(pop, " << convertedInequalityConstraint.ToString(/* oneIndexed = */true) << ");" << std::endl;
        }
        for (auto convertedEqualityConstraint : convertedEqualityConstraints) {
            os << "push!(pop, " << convertedEqualityConstraint.ToString(/* oneIndexed = */true) << ");" << std::endl;
        }
        os << "order=" << hierarchyDegree << ";" << std::endl;
        os << "@CPUtime opt,sol,data=cs_tssos_first(pop,x,order,numeq=" << convertedEqualityConstraints.size() << ",TS=\"MD\",MomentOne=true,solution=true);" << std::endl;

        if (&os != &std::cout) {
            delete(&os);
        }
    }

    void ProblemData::OutputGMSFile(std::string fileName)
    {
        std::ostream &os = (fileName != "" ? *(new std::ofstream(fileName)) : std::cout);

        int n = objectiveFunction.maxIndex, n2 = n + 1;
        for (auto &convertedIndexSet : convertedIndexSets) {
            n = std::max(n, *std::max_element(convertedIndexSet.begin(), convertedIndexSet.end()));
        }
        n += 1;

        os << "Variables objvar,";
        for (int i = 1; i <= n; i++) {
            os << "x" << i << ",;"[i == n];
        }
        os << "\n";

        if (addVariableNonnegativity) {
            os << "Positive Variables ";
            for (int i = 1; i <= n2; i++) {
                os << "x" << i << ",;"[i == n2];
            }
            os << "\n";
        }

        int m = 1 + convertedEqualityConstraints.size() + convertedInequalityConstraints.size();
        os << "Equations ";
        for (int i = 1; i <= m; i++) {
            os << "e" << i << ",;"[i == m];
        }
        os << "\n";
        
        Polynomial neg = -objectiveFunction;
        os << "e1.. objvar" << neg.ToStringGMS(/* withSign = */true) << " =E= 0;" << std::endl;

        int ctr = 1;
        for (int i = 0; i < convertedEqualityConstraints.size(); i++) {
            os << "e" << ++ctr << ".. " << convertedEqualityConstraints[i].ToStringGMS() << " =E= 0;" << std::endl;
        }

        for (int i = 0; i < convertedInequalityConstraints.size(); i++) {
            os << "e" << ++ctr << ".. " << convertedInequalityConstraints[i].ToStringGMS() << " =G= 0;" << std::endl;
        }

        if (enableOriginalVariableNormBound) {
            for (int i = 1; i <= n; i++) {
                if (!addVariableNonnegativity || i >= n2 + 1) {
                    os << "x" << i << ".lo = " << -originalVariableNormBound << ";" << std::endl;
                }
                os << "x" << i << ".up = " << originalVariableNormBound << ";" << std::endl;
            }
        }

        if (&os != &std::cout) {
            delete(&os);
        }

        std::ostream &clique_os = (fileName != "" ? *(new std::ofstream(regex_replace(fileName, std::regex("\\.gms"), "-clique.txt"))) : std::cout);
        clique_os << convertedIndexSets.size() << std::endl;
        for (auto convertedIndexSet : convertedIndexSets) {
            clique_os << convertedIndexSet.size() << std::endl;
            int cnt = 0;
            for (auto elem : convertedIndexSet) {
                if (cnt > 0) {
                    clique_os << " ";
                }
                clique_os << elem + 1;
                cnt++;
            }
            clique_os << std::endl;
        }

        if (&clique_os != &std::cout) {
            delete(&clique_os);
        }
    }

    void ProblemData::OutputMatFile(std::string fileName)
    {
        #ifdef __BUILD_WITH_MATLAB__
        MATFile *pmat = matOpen(fileName.c_str(), "w");
        if (pmat == NULL) {
            std::cerr << "Error creating file " << fileName << std::endl;
            exit(EXIT_FAILURE);
        }

        mxArray *matrixA = mxCreateSparse(A.rows(), A.cols(), A.nonZeros(), mxREAL);
        mwIndex *ir, *jc;
        double *pr;
        ir = mxGetIr(matrixA);
        jc = mxGetJc(matrixA);
        pr = mxGetDoubles(matrixA);

        std::vector<int> AColPerm(A.cols());
        AColPerm[0] = 0;
        int ctr = 1;
        int IndOffset = 1;
        int freeConeSize = 1;
        for (int i = 0; i < psdMatrixSizes.size(); i++) {
            IndOffset += psdMatrixSizes[i] * (psdMatrixSizes[i] + 1) / 2;
        }
        for (int i = 0; i < symmetricMatrixSizes.size(); i++) {
            int symSize = symmetricMatrixSizes[i] * (symmetricMatrixSizes[i] + 1) / 2;
            for (int j = 0; j < symSize; j++) {
                AColPerm[ctr] = IndOffset;
                IndOffset++; ctr++;
            }
            freeConeSize += symSize;
        }
        IndOffset = 1;
        for (int i = 0; i < psdMatrixSizes.size(); i++) {
            int psdSize = psdMatrixSizes[i] * (psdMatrixSizes[i] + 1) / 2;
            for (int j = 0; j < psdSize; j++) {
                AColPerm[ctr] = IndOffset;
                IndOffset++; ctr++;
            }
        }

        Eigen::SparseMatrix<double> A2 = A;
        int ind2 = 0;
        for (int i = 0; i < A2.outerSize(); i++) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A2, i); it; ++it) {
                while (ind2 < sq2Cols.size() && sq2Cols[ind2] < it.col()) ind2++;
                double scalingCoef = (ind2 < sq2Cols.size() && sq2Cols[ind2] == it.col()) ? sqrt(0.5) : 1;
                it.valueRef() = it.value() / (scalingCoef * D[it.row()] * E[it.col()] * scalingFactor);
            }
        }
        
        std::vector<Eigen::Triplet<double>> tripletList;
        int nzCnt = 0;
        for (int i = 0; i < A2.outerSize(); i++) {
            jc[i] = nzCnt;
            for (Eigen::SparseMatrix<double>::InnerIterator it(A2, AColPerm[i]); it; ++it) {
                pr[nzCnt] = it.value();
                ir[nzCnt] = it.row();
                nzCnt++;
            }
        }
        jc[A.outerSize()] = nzCnt;
        
        if (matPutVariable(pmat, "A", matrixA) != 0) {
            std::cerr << __FILE__ << " : " << "Error using matPutvariable on line " << __LINE__ << std::endl;
            exit(EXIT_FAILURE);
        }
        mxDestroyArray(matrixA);

        mxArray *vectorB = mxCreateDoubleMatrix(b.size(), 1, mxREAL);
        for (int i = 0; i < b.size(); i++) {
             mxGetDoubles(vectorB)[i] = (double)b[i] / (D[i] * dualScaler * scalingFactor);
        }
        if (matPutVariable(pmat, "b", vectorB) != 0) {
            std::cerr << __FILE__ << " : " << "Error using matPutvariable on line " << __LINE__ << std::endl;
            exit(EXIT_FAILURE);
        }
        mxDestroyArray(vectorB);

        mxArray *vectorC = mxCreateSparse(c.rows(), c.cols(), c.nonZeros(), mxREAL);
        ir = mxGetIr(vectorC);
        jc = mxGetJc(vectorC);
        pr = mxGetDoubles(vectorC);
        nzCnt = 0;
        for (int i = 0; i < c.outerSize(); i++) {
            jc[i] = nzCnt;
            for (Eigen::SparseVector<double>::InnerIterator it(c, i); it; ++it) {
                pr[nzCnt] = it.value() / (E[it.row()] * primalScaler * scalingFactor);
                ir[nzCnt] = it.row();
                nzCnt++;
            }
        }
        jc[c.outerSize()] = nzCnt;
        if (matPutVariable(pmat, "c", vectorC) != 0) {
            std::cerr << __FILE__ << " : " << "Error using matPutvariable on line " << __LINE__ << std::endl;
            exit(EXIT_FAILURE);
        }
        mxDestroyArray(vectorC);
        
        const char *fieldNames[] = {"f", "s"};
        mxArray *coneInfo = mxCreateStructMatrix(1, 1, 2, fieldNames);
        
        int freeConeID = mxGetFieldNumber(coneInfo, "f");
        mxArray *freeConeFieldValue = mxCreateDoubleMatrix(1, 1, mxREAL);
        *mxGetDoubles(freeConeFieldValue) = freeConeSize * 1.0;
        mxSetFieldByNumber(coneInfo, 0, freeConeID, freeConeFieldValue);
        
        int psdConeID = mxGetFieldNumber(coneInfo, "s");
        mxArray *psdConeFieldValue = mxCreateDoubleMatrix(1, psdMatrixSizes.size(), mxREAL);
        for (int i = 0; i < psdMatrixSizes.size(); i++) {
            mxGetDoubles(psdConeFieldValue)[i] = psdMatrixSizes[i] * 1.0;
        }
        mxSetFieldByNumber(coneInfo, 0, psdConeID, psdConeFieldValue);
        
        if (matPutVariable(pmat, "K", coneInfo) != 0) {
            std::cerr << __FILE__ << " : " << "Error using matPutvariable on line " << __LINE__ << std::endl;
            exit(EXIT_FAILURE);
        }
        mxDestroyArray(coneInfo);

        if (matClose(pmat) != 0) {
            std::cerr << "Error closing file " << fileName << std::endl;
            exit(EXIT_FAILURE);
        }
        #endif //__BUILD_WITH_MATLAB__
    }

    void ProblemData::_GenerateSolutions(int cur, std::vector<std::vector<std::vector<double>>> &solutionCandidates, std::vector<std::vector<double>> &tmp, std::vector<std::vector<double>> &answers)
    {
        if (cur == solutionCandidates.size()) {
	        std::vector<double> ansCand(tmp.size());
	        for (int i = 0; i < ansCand.size(); i++) {
		        if (tmp[i].size() == 0) return;
		        ansCand[i] = std::accumulate(tmp[i].begin(), tmp[i].end(), 0.0) * 1. / tmp[i].size();
	        }
            answers.emplace_back(ansCand);
            return;
        }
	
	    std::vector<int> addInd;

        for (auto solutionCandidate : solutionCandidates[cur]) {
            bool compatible = true;
            for (int i = 0; i < solutionCandidate.size(); i++) {
                int ind = convertedIndexSets[cur][i];
                if (tmp[ind].size() && std::abs(tmp[ind][0] - solutionCandidate[i]) / std::max({1.0, std::abs(tmp[ind][0]), std::abs(solutionCandidate[i])}) >= 1e-3) {
                    compatible = false;
                    break;
                }
                else {
                    tmp[ind].emplace_back(solutionCandidate[i]);
		            addInd.emplace_back(ind);
                }
            }
            if (compatible) {
                _GenerateSolutions(cur + 1, solutionCandidates, tmp, answers);
            }
            for (int i = 0; i < addInd.size(); i++) {
		        tmp[addInd[i]].pop_back();
	        }
      	    addInd.clear();
        }
    }

    std::vector<std::vector<double>> ProblemData::ExtractSolutionsFrom(std::vector<double> &v, bool tellObjVal, double objVal)
    {
        int n = objectiveFunction.maxIndex + (gradientConstraintType != 2);

        std::vector<std::vector<double>> ret;

        using Point = std::vector<double>;
        std::vector<std::vector<Point>> solutionCandidates(convertedIndexSets.size());

        if (addFirstOrderFullMomentMatrix) {
            Eigen::MatrixXd MomentMatrix = Eigen::MatrixXd::Zero(n + 1, n + 1);
            for (int i = 0; i <= n; i++) {
                for (int j = i; j <= n; j++) {
                    Term t = (i == 0 ? (j == 0 ? Term({}) : Term({j - 1})) : Term({i - 1, j - 1}));
                    MomentMatrix(i, j) = MomentMatrix(j, i) = v[termToInteger[t]];
                }
            }
            if (verbose) {
                //std::cout << "Full-order moment matrix of degree 1 : " << std::endl;
                //std::cout << MomentMatrix << std::endl;
            }
            Eigen::ColPivHouseholderQR<Eigen::MatrixXd> QRDecomp(MomentMatrix);
            QRDecomp.setThreshold(1e-5);
            if (verbose) std::cout << "rank(M_{n,1}(y)) = " << QRDecomp.rank() << std::endl;
            if (QRDecomp.rank() == 1) {
                std::vector<double> res(n, 0);
                for (int i = 1; i <= n; i++) {
                    res[i - 1] = MomentMatrix(0, i);
                }
                ret.emplace_back(res);
                return ret;
            }
            if (verbose) std::cout << "rank(M_{n,1})(y) > 1, so we'll continue the ordinary procedure." << std::endl;
        }

        for (int i = 0; i < convertedIndexSets.size(); i++) {
            int originalVariableNum = std::lower_bound(convertedIndexSets[i].begin(), convertedIndexSets[i].end(), n) - convertedIndexSets[i].begin();
            int sz = 0;
            std::vector<std::pair<int, int>> newOrd;
            std::vector<bool> isOriginalTerms(convertedTermSets[i].size(), false);
            std::vector<int> inds;

            for (int j = 0; j < convertedTermSets[i].size(); j++) {
                bool orig = true;
                for (int k = 0; k < convertedTermSets[i][j].size(); k++) orig &= (convertedTermSets[i][j][k] < n);
                isOriginalTerms[j] = orig;
                if (orig) {
                    sz++;
                    newOrd.emplace_back(std::make_pair(convertedTermSets[i][j].size(), j));
                }
            }

            std::sort(newOrd.begin(), newOrd.end());

            for (int j = 0; j < newOrd.size(); j++) {
                inds.emplace_back(newOrd[j].second);
            }

            Eigen::MatrixXd MomentMatrix = Eigen::MatrixXd::Zero(sz, sz);
            for (int j = 0; j < sz; j++) {
                for (int k = j; k < sz; k++) {
                    Term tMerged;
                    Term t1 = convertedTermSets[i][inds[j]], t2 = convertedTermSets[i][inds[k]];
                    std::merge(t1.begin(), t1.end(), t2.begin(), t2.end(), std::back_inserter(tMerged));
                    MomentMatrix(j, k) = v[termToInteger[tMerged]];
                    if (j != k) {
                        MomentMatrix(k, j) = MomentMatrix(j, k);
                    }
                }
            }

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(MomentMatrix);
            if (eigenSolver.info() != Eigen::Success) {
                if (verbose) std::cout << "Error occured on eigen decomposition :(" << std::endl;
                exit(1);
            }
            if (verbose) {
                //std::cout << "eigenvalues : " << std::endl;
                //std::cout << eigenSolver.eigenvalues() << std::endl;
            }
            int positiveNum = 0;
            for (int i = sz - 1; i >= 0; i--) {
                double lambda = eigenSolver.eigenvalues()(i);
                if (lambda / eigenSolver.eigenvalues()(sz - 1) >= 1e-2) {
                    positiveNum++;
                }
                else {
                    break;
                }
            }
            
            Eigen::MatrixXd V = Eigen::MatrixXd::Zero(sz, positiveNum);
            {
                int ptr = 0;
                for (int i = 0; i < positiveNum; i++) {
                    double lambda = eigenSolver.eigenvalues()(sz - i - 1);
                    V.col(ptr) = sqrt(lambda) * eigenSolver.eigenvectors().col(sz - i - 1);
                    ptr++;
                }
            }

            int colNumOfV = V.cols();
            std::vector<Term> tmpTerms;

            for (int j = 0; j < newOrd.size(); j++) {
                tmpTerms.emplace_back(convertedTermSets[i][inds[j]]);
            }

            std::vector<Term> baseTerms;
            {
                int ptr = 0;
                for (int j = 0; j < colNumOfV; ) {
                    if (ptr == newOrd.size()) {
                        if (verbose) std::cout << "[ExtractSolution] [Error] making column echelon form failed" << std::endl;
                        return std::vector<std::vector<double>>();
                    }

                    int bigCol = -1;
                    for (int k = j; k < colNumOfV; k++) {
                        if (bigCol == -1 || abs(V(ptr, k)) > abs(V(ptr, bigCol))) {
                            bigCol = k;
                        }
                    }

                    if (abs(V(ptr, bigCol)) < 1e-5) {
                        for (int k = j; k < colNumOfV; k++) V(ptr, k) = 0;
                        ptr++; continue;
                    }
                    else if (bigCol != j) {
                        V.col(j).swap(V.col(bigCol));
                    }

                    double rowMax = abs(V(ptr, j));
                    for (int k = ptr + 1; k < newOrd.size(); k++) rowMax = std::max(rowMax, abs(V(k, j)));
                    while (ptr < newOrd.size() && abs(V(ptr, j)) < 1e-5 * rowMax) ptr++;

                    if (ptr == newOrd.size() || tmpTerms[ptr].size() == hierarchyDegree) {
                        if (verbose) std::cout << "[ExtractSolution] [Error] degree of the monomial is too large" << std::endl;
                        if (verbose) std::cout << V << std::endl;
                        return std::vector<std::vector<double>>();
                    }

                    baseTerms.emplace_back(tmpTerms[ptr]);

                    V.col(j) /= V(ptr, j);
                    for (int k = 0; k < colNumOfV; k++) {
                        if (j != k) {
                            V.col(k) -= V(ptr, k) * V.col(j);
                            V(ptr, k) = 0;
                        }
                    }
                    j++; ptr++;
                }
            }

            std::map<Term, int> t2i;
            for (int j = 0; j < sz; j++) {
                t2i[tmpTerms[j]] = j;
            }

            //make a Multiplication Matrix
            std::random_device seed_gen;
            std::default_random_engine engine(seed_gen());

            std::uniform_real_distribution<> dist(0.0, 1.0);
            std::vector<double> lambdas(originalVariableNum);
            for (int j = 0; j < originalVariableNum; j++) {
                lambdas[j] = 1 - dist(engine);
            }
            double tot = std::accumulate(lambdas.begin(), lambdas.end(), 0.0);
            for (int j = 0; j < originalVariableNum; j++) {
                lambdas[j] /= tot;
            }

            Eigen::MatrixXd N = Eigen::MatrixXd::Zero(colNumOfV, colNumOfV);
            std::vector<Eigen::MatrixXd> Ns;
            for (int j = 0; j < originalVariableNum; j++) {
                Eigen::MatrixXd Nj(colNumOfV, colNumOfV);
                for (int k = 0; k < colNumOfV; k++) {
                    Term term = baseTerms[k]; term.emplace_back(convertedIndexSets[i][j]);
                    std::sort(term.begin(), term.end());
                    Nj.row(k) = V.row(t2i[term]);
                }
                N += lambdas[j] * Nj;
                Ns.emplace_back(Nj);
            }
            
            Eigen::RealSchur<Eigen::MatrixXd> NSchur(N);
            if (NSchur.matrixT().isUpperTriangular(1e-5) == false) {
                if (verbose) std::cout << "[ExtractSolution] [Error] Schur Decompostion contains a complex number (i = " << i << ")" << std::endl;
                if (verbose) std::cout << NSchur.matrixT() << std::endl;
                return std::vector<std::vector<double>>();
            }

            solutionCandidates[i].resize(colNumOfV);
            for (int j = 0; j < colNumOfV; j++) solutionCandidates[i][j].resize(originalVariableNum);

            Eigen::MatrixXd Q = NSchur.matrixU();
            for (int j = 0; j < originalVariableNum; j++) {
                for (int k = 0; k < colNumOfV; k++) {
                    solutionCandidates[i][k][j] = Q.col(k).transpose() * Ns[j] * Q.col(k);
                }
            }
            if (verbose) {
                for (int k = 0; k < colNumOfV; k++) {
                    for (int j = 0; j < originalVariableNum; j++) {
                        std::cout << convertedIndexSets[i][j] << " = " << solutionCandidates[i][k][j] << ",";
                    }
                    std::cout << std::endl;
                }
            }
        }

        std::vector<std::vector<double>> tmp(n);

        _GenerateSolutions(0, solutionCandidates, tmp, ret);
	
	    int cnt = 1;
        for (auto it = ret.begin(); it != ret.end(); ) {
            if (verbose) {
                std::cout << "candidate " << cnt << std::endl; cnt++;
                Analyze(*it);
            }
            std::vector<double> x = *it;
            double fx = objectiveFunction.Evaluate(x);

            if (tellObjVal == true && std::abs(fx - objVal) >= 1e-5) {
                it = ret.erase(it); continue;
            }

            double norm = 0;
            for (int i = 0; i < x.size(); i++) norm += x[i] * x[i];
            if (std::abs(norm - 1) >= 1e-5) {
                it = ret.erase(it); continue;
            }

            double nablaFNorm = 0;
            std::vector<double> dx(x.size());
            for (int i = 0; i < x.size(); i++) {
                dx[i] = objectiveFunction.DifferentiateBy(i).Evaluate(x);
                nablaFNorm += dx[i] * dx[i];
            }

            double innerProd = 0;
            for (int i = 0; i < x.size(); i++) innerProd += x[i] * dx[i];
            if (std::abs(std::abs(innerProd) / (sqrt(norm) * sqrt(nablaFNorm)) - 1) >= 1e-5) {
                it = ret.erase(it); continue;
            }
            it++;
        }

        return ret;
    }

    double ProblemData::Apply(std::vector<double> &x)
    {
        return objectiveFunction.Evaluate(x);
    }

    void ProblemData::Analyze(std::vector<double> &x)
    {
        std::cout << std::resetiosflags(std::ios_base::floatfield);
        std::cout << "(";
        for (int i = 0; i < x.size(); i++) {
            std::cout << x[i];
            if (i + 1 == x.size()) {
                std::cout << ")" << std::endl;
            }
            else {
                std::cout << ", ";
            }
        }

        std::cout << "f(x) = " << objectiveFunction.Evaluate(x) << std::endl;

        std::cout << "||x|| = ";
        double norm = 0;
        for (int i = 0; i < x.size(); i++) norm += x[i] * x[i];
        std::cout << sqrt(norm) << std::endl;


        double nablaFNorm = 0;

        for (int i = 0; i < x.size(); i++) {
            double dx_i = objectiveFunction.DifferentiateBy(i).Evaluate(x);
            nablaFNorm += dx_i * dx_i;
        }
        std::cout << "||\\nabla f|| = " << sqrt(nablaFNorm) << std::endl;

        double innerProd = 0;
        for (int i = 0; i < x.size(); i++) innerProd += x[i] * objectiveFunction.DifferentiateBy(i).Evaluate(x);
        std::cout << "cos(\\nabla f, x) = " << innerProd / (sqrt(norm) * sqrt(nablaFNorm)) << std::endl;
    }
}
