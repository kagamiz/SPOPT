#include "Solver.hpp"

namespace SPOPT {

    void Solver::LoadConfig(std::string fileName)
    {
        solverConfig = YAML::LoadFile(fileName);
        basicParam.enableAndersonAcceleration = solverConfig["enableAndersonAcceleration"].as<bool>(false);
        basicParam.primalTolerance            = solverConfig["primalTolerance"].as<double>(1e-3);
        basicParam.dualTolerance              = solverConfig["dualTolerance"].as<double>(1e-3);
        basicParam.gapTolerance               = solverConfig["gapTolerance"].as<double>(1e-3);
        basicParam.reportFrequency            = solverConfig["reportFrequency"].as<int>(1);
    }

    void Solver::Solve(const ProblemData &problemData)
    {
        SetUpFrom(problemData);
        int iterationCounter = 0;
        Eigen::VectorXd v_i = ConstructInitialPoint(problemData);

        if (basicParam.reportFrequency > 0) {
            std::cout << "Ite No.  PFEAS   DFEAS    GAP         PVAL            DVAL     " << std::endl;
            std::cout << "===============================================================" << std::endl;
        }

        while (!IsTerminationCriterionSatisfied(problemData, v_i)) {
            if (basicParam.reportFrequency > 0 && iterationCounter % basicParam.reportFrequency == 0) {
                std::cout << std::setw(7)  << iterationCounter << " "
                          << std::setw(7)  << std::scientific << GetPrimalInfeasibility(problemData, v_i) << " "
                          << std::setw(7)  << std::scientific << GetDualInfeasibility(problemData, v_i)   << " "
                          << std::setw(7)  << std::scientific << GetGap(problemData, v_i)                 << " "
                          << std::setw(15) << std::scientific << GetPrimalObjValue(problemData, v_i)      << " "
                          << std::setw(15) << std::scientific << GetDualObjValue(problemData, v_i)
                          << std::endl;
                std::cout << std::resetiosflags;
            }
            iterationCounter++;
            v_i = ApplyFixedPointFunction(problemData, v_i);
            UpdateParameter(problemData, v_i);
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