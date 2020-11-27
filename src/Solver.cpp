#include "Solver.hpp"

namespace SPOPT {

    Solver::Solver(std::string fileName)
    {
        LoadConfig(fileName);
    }

    void Solver::LoadConfig(std::string fileName)
    {
        solverConfig = YAML::LoadFile(fileName);

        param.enableAndersonAcceleration = solverConfig["enableAndersonAcceleration"].as<bool>(false);
        param.primalTolerance            = solverConfig["primalTolerance"].as<double>(1e-3);
        param.dualTolerance              = solverConfig["dualTolerance"].as<double>(1e-3);
        param.gapTolerance               = solverConfig["gapTolerance"].as<double>(1e-3);
    }

    void Solver::Solve(const ProblemData &problemData)
    {
        int iterationCounter = 0;
        Eigen::VectorXd v_i = ConstructInitialPoint(problemData);

        while (!IsTerminationCriterionSatisfied(problemData, v_i)) {
            iterationCounter++;
            v_i = ApplyFixedPointFunction(problemData, v_i);
        }
    }
}