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
    }

    void Solver::ShowHeader()
    {
        std::cout << "Ite No.   PFEAS    DFEAS     GAP         PVAL             DVAL      " << std::endl;
        std::cout << "====================================================================" << std::endl;
    }

    void Solver::ShowIterInfo(int iterID, const ProblemData &problemData, Eigen::VectorXd &v)
    {
        std::cout << std::setw(7)  << iterID << " "
                  << std::setw(7)  << std::setprecision(2)  << std::scientific << GetPrimalInfeasibility(problemData, v) << " "
                  << std::setw(7)  << std::setprecision(2)  << std::scientific << GetDualInfeasibility(problemData, v)   << " "
                  << std::setw(7)  << std::setprecision(2)  << std::scientific << GetGap(problemData, v)                 << " "
                  << std::setw(15) << std::setprecision(10) << std::scientific << GetPrimalObjValue(problemData, v)      << " "
                  << std::setw(15) << std::setprecision(10) << std::scientific << GetDualObjValue(problemData, v)
                  << std::endl;
        std::cout << std::resetiosflags(std::ios_base::floatfield);
    }

    void Solver::Solve(const ProblemData &problemData)
    {
        SetUpFrom(problemData);
        int iterationCounter = 0;
        Eigen::VectorXd v_i = ConstructInitialPoint(problemData);

        AaWork *aawork;
        if (basicParam.enableAndersonAcceleration) {
            aawork = aa_init(v_i.size(), basicParam.AAMemoryLength, basicParam.AAType, basicParam.AAEta);
        }

        if (basicParam.reportFrequency > 0) {
            ShowHeader();
        }
        while (!IsTerminationCriterionSatisfied(problemData, v_i)) {
            if (basicParam.reportFrequency > 0 && iterationCounter % basicParam.reportFrequency == 0) {
                ShowIterInfo(iterationCounter, problemData, v_i);
            }
            iterationCounter++;
            Eigen::VectorXd v_prv;
            if (basicParam.enableAndersonAcceleration) {
                v_prv = v_i;
            }
            v_i = ApplyFixedPointFunction(problemData, v_i);
            if (basicParam.enableAndersonAcceleration) {
                double cur[v_i.size()], prv[v_i.size()];
                Eigen::VectorXd::Map(cur, v_i.rows()) = v_i;
                Eigen::VectorXd::Map(prv, v_prv.rows()) = v_prv;
                aa_apply(cur, prv, aawork);
                for (int i = 0; i < v_i.size(); i++) {
                    v_i(i) = cur[i];
                }
            }
            UpdateParameter(problemData, v_i);
        }
        ShowIterInfo(iterationCounter, problemData, v_i);

        std::cout << "Solver terminated on iteration #" << iterationCounter << "!!" << std::endl;
    }

    bool Solver::IsTerminationCriterionSatisfied(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        double primalInfeasibility = GetPrimalInfeasibility(problemData, v);
        double dualInfeasibility   = GetDualInfeasibility(problemData, v);
        double gap                 = GetGap(problemData, v);

        return primalInfeasibility <= basicParam.primalTolerance && dualInfeasibility <= basicParam.dualTolerance && gap <= basicParam.gapTolerance;
    }
}