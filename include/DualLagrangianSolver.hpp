#ifndef __DUAL_LAGRANGIAN_SOLVER_HPP__
#define __DUAL_LAGRANGIAN_SOLVER_HPP__

#include "Solver.hpp"

namespace SPOPT {

    class DualLagrangianParameter {
        public:
            const double muMin = 1e-10;
            const double muMax = 1e+10;
            double mu;
            double rho;
            double gamma;
            int h1, h2, h3, h4;
    };

    class DualLagrangianSolver : public Solver {
        public:
            // constructors
            DualLagrangianSolver(){};
            // Reads solver settings from the file `filename` and store it on `solverConfig`
            // when instanciating the Solver.
            // Internally, it calls a function LoadConfig.
            DualLagrangianSolver(std::string fileName);

            void LoadConfig(std::string fileName);

        private:

            void SetUpFrom(const ProblemData &problemData);
            Eigen::VectorXd ConstructInitialPoint(const ProblemData &problemData);
            Eigen::VectorXd ApplyFixedPointFunction(const ProblemData &problemData, const Eigen::VectorXd &v);
            void UpdateParameter(const ProblemData &problemData, const Eigen::VectorXd &v);

            Eigen::VectorXd CalcInvAAt(const ProblemData &problemData, Eigen::VectorXd &v);

            double GetPrimalInfeasibility(const ProblemData &problemData, const Eigen::VectorXd &v);
            double GetDualInfeasibility(const ProblemData &problemData, const Eigen::VectorXd &v);
            double GetPrimalObjValue(const ProblemData &problemData, const Eigen::VectorXd &v) = 0;
            double GetDualObjValue(const ProblemData &problemData, const Eigen::VectorXd &v) = 0;
            double GetGap(const ProblemData &problemData, const Eigen::VectorXd &v);

            DualLagrangianParameter dualParam;

            Eigen::SparseMatrix<double> At;
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> AAt;
            Eigen::VectorXd sparseAAt;

            int it_pinf, it_dinf;
    };
}
#endif //__DUAL_LAGRANGIAN_SOLVER_HPP__