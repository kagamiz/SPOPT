#include "Solver.hpp"

namespace SPOPT {

    void Solver::LoadConfig(std::string fileName)
    {
        solverConfig = YAML::LoadFile(fileName);
        basicParam.enableAndersonAcceleration = solverConfig["enableAndersonAcceleration"].as<bool>(false);
        basicParam.AAType                     = solverConfig["AAType"].as<bool>(true);
        basicParam.AAMemoryLength             = solverConfig["AAMemoryLength"].as<int>(10);
        basicParam.AAEta                      = solverConfig["AAEta"].as<double>(1e-8);
        basicParam.primalTolerance            = solverConfig["primalTolerance"].as<double>(1e-3);
        basicParam.dualTolerance              = solverConfig["dualTolerance"].as<double>(1e-3);
        basicParam.gapTolerance               = solverConfig["gapTolerance"].as<double>(1e-3);
        basicParam.reportFrequency            = solverConfig["reportFrequency"].as<int>(1);
        basicParam.iterationLimit             = solverConfig["iterationLimit"].as<int>(2500);
    }

    void Solver::ShowHeader()
    {
        std::cout << "Ite No.   PFEAS    DFEAS     GAP         PVAL             DVAL        TIME      " << std::endl;
        std::cout << "================================================================================" << std::endl;
    }

    void Solver::ShowIterInfo(int iterID, const ProblemData &problemData, Eigen::VectorXd &v)
    {
        std::cout << std::setw(7)  << iterID << " "
                  << std::setw(7)  << std::setprecision(2)  << std::scientific << GetPrimalInfeasibility(problemData, v) << " "
                  << std::setw(7)  << std::setprecision(2)  << std::scientific << GetDualInfeasibility(problemData, v)   << " "
                  << std::setw(7)  << std::setprecision(2)  << std::scientific << GetGap(problemData, v)                 << " "
                  << std::setw(15) << std::setprecision(10) << std::scientific << GetPrimalObjValue(problemData, v)      << " "
                  << std::setw(15) << std::setprecision(10) << std::scientific << GetDualObjValue(problemData, v)        << " "
                  << std::setw(15) << std::setprecision(10) << std::scientific << elapsedTime << " "
                  << std::endl;
        std::cout << std::resetiosflags(std::ios_base::floatfield);
    }

    std::pair<double, std::vector<double>> Solver::Solve(const ProblemData &problemData)
    {
        SetUpFrom(problemData);
        int iterationCounter = 0;
        Eigen::VectorXd v_i = ConstructInitialPoint(problemData);
        Eigen::VectorXd v_prv;

        AaWork *aawork;
        double *aacur, *aaprv;
        if (basicParam.enableAndersonAcceleration) {
            aawork = aa_init(v_i.size(), basicParam.AAMemoryLength, basicParam.AAType, basicParam.AAEta);
            aacur = (double *)malloc(v_i.size() * sizeof(double));
            aaprv = (double *)malloc(v_i.size() * sizeof(double));
        }
        if (basicParam.reportFrequency > 0) {
            ShowHeader();
        }
        while (iterationCounter < basicParam.iterationLimit && !IsTerminationCriterionSatisfied(problemData, v_i)) {
            clock_t start = clock();
            if (iterationCounter > 0 && basicParam.enableAndersonAcceleration) {
                Eigen::VectorXd::Map(aacur, v_i.size())   = v_i;
                Eigen::VectorXd::Map(aaprv, v_prv.size()) = v_prv;
                aa_apply(aacur, aaprv, aawork);
                for (int i = 0; i < v_i.size(); i++) {
                    v_i(i) = aacur[i];
                }
            }

            if (basicParam.enableAndersonAcceleration) {
                v_prv = v_i;
            }
            
            v_i = ApplyFixedPointFunction(problemData, v_i);
            
            UpdateParameter(problemData, v_i);
            clock_t end = clock();
            elapsedTime += (double)(end - start) / CLOCKS_PER_SEC;

            if (basicParam.reportFrequency > 0 && iterationCounter % basicParam.reportFrequency == 0) {
                ShowIterInfo(iterationCounter, problemData, v_i);
            }
            iterationCounter++;
        }
        ShowIterInfo(iterationCounter, problemData, v_i);

        std::cout << "Solver terminated on iteration #" << iterationCounter << "!!" << std::endl;

        if (basicParam.enableAndersonAcceleration) {
            aa_finish(aawork);
            free(aacur);
            free(aaprv);
        }

        return std::make_pair(0, std::vector<double>());
    }

    void Solver::EnumerateStationaryPoints(std::string fileName)
    {
        clock_t start = clock();
        YAML::Node problemDataConfig = YAML::LoadFile(fileName);
        
        problemDataConfig["enableUpperBound"] = false;
        problemDataConfig["enableLowerBound"] = false;
        problemDataConfig["gradientConstraintType"] = 0;

        int originalDegree = problemDataConfig["hierarchyDegree"].as<int>();

        double lambda_min;
        std::vector<std::vector<double>> min_vecs;
        while (true) {
            ProblemData p;
            p.LoadConfig(problemDataConfig);
            p.ConstructSDP();
            auto [min_cand, tms] = Solve(p);
            min_vecs = p.ExtractSolutionsFrom(tms, true, min_cand);

            if (min_vecs.size()) {
                lambda_min = min_cand;
                if (problemDataConfig["verbose"].as<bool>(false)) std::cout << "minimum eigvalue extracted! hierarchyDegree = " << problemDataConfig["hierarchyDegree"].as<int>() << std::endl;
                break;
            }
            else {
                problemDataConfig["hierarchyDegree"] = problemDataConfig["hierarchyDegree"].as<int>() + 1;
                if (problemDataConfig["hierarchyDegree"].as<int>() >= 6) {
                    std::cout << "[fatal] Stopping Extraction Process because relaxation stage is too large." << std::endl;
                    return;
                }
            }
        }


        problemDataConfig["hierarchyDegree"] = originalDegree;
        problemDataConfig["maximize"] = true;
        double lambda_max;
        std::vector<std::vector<double>> max_vecs;
        while (true) {
            ProblemData p;
            p.LoadConfig(problemDataConfig);
            p.ConstructSDP();
            auto [max_cand, tms] = Solve(p);
            max_vecs = p.ExtractSolutionsFrom(tms, true, max_cand);

            if (max_vecs.size()) {
                lambda_max = -max_cand;
                if (problemDataConfig["verbose"].as<bool>(false)) std::cout << "maximum eigvalue extracted! hierarchyDegree = " << problemDataConfig["hierarchyDegree"].as<int>() << std::endl;
                break;
            }
            else {
                problemDataConfig["hierarchyDegree"] = problemDataConfig["hierarchyDegree"].as<int>() + 1;
                if (problemDataConfig["hierarchyDegree"].as<int>() >= 6) {
                    std::cout << "[fatal] Stopping Extraction Process because relaxation stage is too large." << std::endl;
                    return;
                }
            }
        }

        std::map< double, std::vector<std::vector<double>> > realEigenPairs;

        realEigenPairs[lambda_min] = min_vecs;
        realEigenPairs[lambda_max] = max_vecs;
        
        problemDataConfig["gradientConstraintType"] = 2;

        std::stack<std::pair<double, double>> intervals;
        intervals.emplace(lambda_min, lambda_max);
        while (intervals.size()) {
            auto [lo, hi] = intervals.top(); intervals.pop();

            if (problemDataConfig["verbose"]) std::cout << "[" << lo << ", " << hi << "]" << std::endl;

            double mid = (lo + hi) / 2;

            problemDataConfig["hierarchyDegree"] = originalDegree;
            problemDataConfig["maximize"] = true;
            problemDataConfig["enableUpperBound"] = true;
            problemDataConfig["enableLowerBound"] = false;
            problemDataConfig["upperBoundConstant"] = mid;

            double tmp_lo;
            std::vector<std::vector<double>> tmp_lo_vecs;
            while (true) {
                ProblemData p;
                p.LoadConfig(problemDataConfig);
                p.ConstructSDP();
                auto [min_cand, tms] = Solve(p);
                tmp_lo_vecs = p.ExtractSolutionsFrom(tms, true, min_cand);

                if (tmp_lo_vecs.size()) {
                    tmp_lo = min_cand;
                    if (problemDataConfig["verbose"].as<bool>(false)) std::cout << "minimum eigvalue extracted! hierarchyDegree = " << problemDataConfig["hierarchyDegree"].as<int>() << std::endl;
                    break;
                }
                else {
                    problemDataConfig["hierarchyDegree"] = problemDataConfig["hierarchyDegree"].as<int>() + 1;
                    if (problemDataConfig["hierarchyDegree"].as<int>() >= 6) {
                        std::cout << "Stopping Extraction Process because relaxation stage is too large." << std::endl;
                        return;
                    }
                }
            }

            if (std::abs(tmp_lo - lo) >= 1e-5) {
                realEigenPairs[tmp_lo] = tmp_lo_vecs;
                intervals.emplace(lo, tmp_lo);
            }

            problemDataConfig["hierarchyDegree"] = originalDegree;
            problemDataConfig["maximize"] = false;
            problemDataConfig["enableUpperBound"] = false;
            problemDataConfig["enableLowerBound"] = true;
            problemDataConfig["lowerBoundConstant"] = mid;

            double tmp_hi;
            std::vector<std::vector<double>> tmp_hi_vecs;
            while (true) {
                ProblemData p;
                p.LoadConfig(problemDataConfig);
                p.ConstructSDP();
                auto [max_cand, tms] = Solve(p);
                tmp_hi_vecs = p.ExtractSolutionsFrom(tms, true, max_cand);

                if (max_vecs.size()) {
                    tmp_hi = -max_cand;
                    if (problemDataConfig["verbose"].as<bool>(false)) std::cout << "maximum eigvalue extracted! hierarchyDegree = " << problemDataConfig["hierarchyDegree"].as<int>() << std::endl;
                    break;
                }
                else {
                    problemDataConfig["hierarchyDegree"] = problemDataConfig["hierarchyDegree"].as<int>() + 1;
                    if (problemDataConfig["hierarchyDegree"].as<int>() >= 6) {
                        std::cout << "Stopping Extraction Process because relaxation stage is too large." << std::endl;
                        return;
                    }
                }
            }
            if (std::abs(hi - tmp_hi) >= 1e-5) {
                if (std::abs(tmp_lo - tmp_hi) >= 1e-5) realEigenPairs[tmp_hi] = tmp_hi_vecs;
                intervals.emplace(tmp_hi, hi);
            }
        }

        clock_t end = clock();
        std::cout << "elapsed = " << (double)(end - start) / CLOCKS_PER_SEC << std::endl;

        for (auto it = realEigenPairs.begin(); it != realEigenPairs.end(); it++) {
            std::cout << it->first << " (num of eigenvector : " << it->second.size() << ")" << std::endl;
        }
    }

    bool Solver::IsTerminationCriterionSatisfied(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        double primalInfeasibility = GetPrimalInfeasibility(problemData, v);
        double dualInfeasibility   = GetDualInfeasibility(problemData, v);
        double gap                 = GetGap(problemData, v);

        return primalInfeasibility <= basicParam.primalTolerance && dualInfeasibility <= basicParam.dualTolerance && gap <= basicParam.gapTolerance;
    }
}