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

    void Solver::Solve(const ProblemData &problemData)
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
    }

    bool Solver::IsTerminationCriterionSatisfied(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        double primalInfeasibility = GetPrimalInfeasibility(problemData, v);
        double dualInfeasibility   = GetDualInfeasibility(problemData, v);
        double gap                 = GetGap(problemData, v);

        return primalInfeasibility <= basicParam.primalTolerance && dualInfeasibility <= basicParam.dualTolerance && gap <= basicParam.gapTolerance;
    }
}