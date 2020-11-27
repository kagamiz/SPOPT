#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#include <yaml-cpp/yaml.h>

#include "LinearAlgebra.hpp"
#include "StandardLibraries.hpp"
#include "ProblemData.hpp"

namespace SPOPT {

    class Parameter {
        public:
            bool enableAndersonAcceleration;
            double primalTolerance;
            double dualTolerance;
            double gapTolerance;
    };

    class Solver {
        public:
            // constructors
            Solver(){};
            // Reads solver settings from the file `filename` and store it on `solverConfig`
            // when instanciating the Solver.
            // Internally, it calls a function LoadConfig.
            Solver(std::string fileName);

            // Reads solver settings from the file `filename` and store it on `solverConfig`.
            // TODO: error handling
            virtual void LoadConfig(std::string fileName);

            void Solve(const ProblemData &problemData);
            void EnumerateStationaryPoints(const ProblemData &problemData);
        private:

            virtual Eigen::VectorXd ConstructInitialPoint(const ProblemData &problemData) = 0;
            virtual bool IsTerminationCriterionSatisfied(const ProblemData &problemData, const Eigen::VectorXd &v) = 0;
            virtual Eigen::VectorXd ApplyFixedPointFunction(const ProblemData &problemData, const Eigen::VectorXd &v) = 0;

            YAML::Node solverConfig;
            Parameter param;
    };
}
#endif //__SOLVER_HPP__