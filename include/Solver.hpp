#ifndef __SOLVER_HPP__
#define __SOLVER_HPP__

#include <yaml-cpp/yaml.h>

#include "LinearAlgebra.hpp"
#include "StandardLibraries.hpp"
#include "ProblemData.hpp"
#include "aa.h"
#include "aa_blas.h"

namespace SPOPT {

    class Parameter {
        public:
            bool enableAndersonAcceleration;
            bool AAType;
            int AAMemoryLength;
            double AAEta;

            double primalTolerance;
            double dualTolerance;
            double gapTolerance;
            int reportFrequency;
    };

    class Solver {
        public:
            // Reads solver settings from the file `filename` and store it on `solverConfig`.
            // TODO: error handling
            virtual void LoadConfig(std::string fileName);

            void Solve(const ProblemData &problemData);
            void EnumerateStationaryPoints(const ProblemData &problemData);
        
        protected:
            bool IsTerminationCriterionSatisfied(const ProblemData &problemData, const Eigen::VectorXd &v);

            void ShowHeader();
            void ShowIterInfo(int iterID, const ProblemData &problemData, Eigen::VectorXd &v);

            virtual void SetUpFrom(const ProblemData &problemData) = 0;
            virtual Eigen::VectorXd ConstructInitialPoint(const ProblemData &problemData) = 0;
            virtual Eigen::VectorXd ApplyFixedPointFunction(const ProblemData &problemData, const Eigen::VectorXd &v) = 0;
            virtual void UpdateParameter(const ProblemData &problemData, const Eigen::VectorXd &v) = 0;

            virtual double GetPrimalInfeasibility(const ProblemData &problemData, const Eigen::VectorXd &v) = 0;
            virtual double GetDualInfeasibility(const ProblemData &problemData, const Eigen::VectorXd &v) = 0;
            virtual double GetPrimalObjValue(const ProblemData &problemData, const Eigen::VectorXd &v) = 0;
            virtual double GetDualObjValue(const ProblemData &problemData, const Eigen::VectorXd &v) = 0;
            virtual double GetGap(const ProblemData &problemData, const Eigen::VectorXd &v) = 0;

            inline const Eigen::SparseMatrix<double> &MatrixA(const ProblemData &problemData) { return problemData.A; }

            inline const Eigen::VectorXd &VectorB(const ProblemData &problemData) { return problemData.b; }
            inline const double originalBNorm(const ProblemData &problemData) { return problemData.originalBNorm; }

            inline const Eigen::SparseVector<double> &VectorC(const ProblemData &problemData) { return problemData.c; }
            inline const double originalCNorm(const ProblemData &problemData) { return problemData.originalCNorm; }
            
            inline const std::vector<int> &psdMatrixSizes(const ProblemData &problemData) { return problemData.psdMatrixSizes; }

            inline const double primalScaler(const ProblemData &problemData) { return problemData.primalScaler; }
            inline const double dualScaler(const ProblemData &problemData) { return problemData.dualScaler; }
            inline const Eigen::VectorXd &VectorD(const ProblemData &problemData) { return problemData.D; }
            inline const Eigen::VectorXd &VectorE(const ProblemData &problemData) { return problemData.E; }

            YAML::Node solverConfig;
            Parameter basicParam;
    };
}
#endif //__SOLVER_HPP__