#include "MOSEKSolver.hpp"

using namespace mosek::fusion;
using namespace monty;

namespace SPOPT {

    MOSEKSolver::MOSEKSolver(std::string fileName)
    {
        LoadConfig(fileName);
    }

    void MOSEKSolver::LoadConfig(std::string fileName)
    {
        Solver::LoadConfig(fileName);
    }

    void MOSEKSolver::Solve(const ProblemData &problemData)
    {
        Model::t M = new Model("SOS_Hierarchy"); auto _M = finally([&]() { M->dispose(); } );
        M->setLogHandler([=](const std::string &msg) { std::cout << msg << std::flush; } );

        std::shared_ptr<ndarray<Variable::t>> variables(new ndarray<Variable::t>(1 + psdMatrixSizes(problemData).size() + symmetricMatrixSizes(problemData).size()));
        (*variables)[0] = M->variable();
        
        {
            int ind = 1;
            for (auto &psdMatrixSize : psdMatrixSizes(problemData)) {
                (*variables)[ind++] = M->variable(Domain::isLinPSD(psdMatrixSize));
            }
            for (auto &symmetricMatrixSize : symmetricMatrixSizes(problemData)) {
                (*variables)[ind++] = M->variable(symmetricMatrixSize * (symmetricMatrixSize + 1) / 2);
            }
        }
        Variable::t x = Var::stack(variables, 0);
        M->objective(ObjectiveSense::Maximize, x->index(0));
        
        const Eigen::SparseMatrix<double> &Aref = MatrixA(problemData);
        std::shared_ptr<ndarray<int>> subi(new ndarray<int>(Aref.nonZeros())), subj(new ndarray<int>(Aref.nonZeros()));
        std::shared_ptr<ndarray<double>> vals(new ndarray<double>(Aref.nonZeros()));

        int ind = 0;
        int ind2 = 0;
        for (int i = 0; i < Aref.outerSize(); i++) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(Aref, i); it; ++it) {
                while (ind2 < sq2Cols(problemData).size() && sq2Cols(problemData)[ind2] < it.col()) ind2++;
                double scalingCoef = (ind2 < sq2Cols(problemData).size() && sq2Cols(problemData)[ind2] == it.col()) ? sqrt(0.5) : 1;
                (*subi)[ind] = it.row();
                (*subj)[ind] = it.col();
                (*vals)[ind] = it.value() / (scalingCoef * VectorD(problemData)[it.row()] * VectorE(problemData)[it.col()] * scalingFactor(problemData));
                ind++;
            }
        }

        Matrix::t A = Matrix::sparse((int)Aref.rows(), (int)Aref.cols(), subi, subj, vals);

        const Eigen::VectorXd &bref = VectorB(problemData);
        std::shared_ptr<ndarray<double>> b(new ndarray<double>(Aref.rows()));
        for (int i = 0; i < bref.size(); i++) {
            (*b)[i] = bref(i) / (VectorD(problemData)[i] * dualScaler(problemData) * scalingFactor(problemData));
        }
        M->constraint(Expr::mul(A, x), Domain::equalsTo(b));
        M->solve();

        std::cout << "Solution : " << (*(x->level()))[0] << std::endl;
    }
}