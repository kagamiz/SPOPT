#ifndef __HSDE_SOLVER_HPP__
#define __HSDE_SOLVER_HPP__

#include "Solver.hpp"

namespace SPOPT {

    class HSDEParameter {
        public:
            double alpha;
    };

    class HSDESolver : public Solver {
        public:
            // constructors
            HSDESolver(){};
            // Reads solver settings from the file `filename` and store it on `solverConfig`
            // when instanciating the Solver.
            // Internally, it calls a function LoadConfig.
            HSDESolver(std::string fileName);

            void LoadConfig(std::string fileName);

        private:

            void SetUpFrom(const ProblemData &problemData);
            Eigen::VectorXd ConstructInitialPoint(const ProblemData &problemData);
            Eigen::VectorXd ApplyFixedPointFunction(const ProblemData &problemData, const Eigen::VectorXd &v);
            void UpdateParameter(const ProblemData &problemData, const Eigen::VectorXd &v);

            Eigen::VectorXd CalcMInv(const ProblemData &problemData, Eigen::VectorXd &v);

            double GetPrimalInfeasibility(const ProblemData &problemData, const Eigen::VectorXd &v);
            double GetDualInfeasibility(const ProblemData &problemData, const Eigen::VectorXd &v);
            double GetPrimalObjValue(const ProblemData &problemData, const Eigen::VectorXd &v);
            double GetDualObjValue(const ProblemData &problemData, const Eigen::VectorXd &v);
            double GetGap(const ProblemData &problemData, const Eigen::VectorXd &v);

            std::vector<double> GetDualVariable(const ProblemData &problemData, const Eigen::VectorXd &v);

            HSDEParameter hsdeParam;

            Eigen::VectorXd zeta;
            Eigen::VectorXd MInvZeta;
            double zetaDotMInvZeta;

            Eigen::SparseMatrix<double> At;
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> IAAt;
            Eigen::VectorXd sparseIAAt;
    };
}
#endif //__HSDE_SOLVER_HPP__