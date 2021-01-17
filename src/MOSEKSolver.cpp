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

    std::pair<double, std::vector<double>> MOSEKSolver::Solve(const ProblemData &problemData)
    {
        Model::t M = new Model("SOS_Hierarchy"); auto _M = finally([&]() { M->dispose(); } );
        if (verbose(problemData)) M->setLogHandler([=](const std::string &msg) { std::cout << msg << std::flush; } );

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
        M->constraint("linEq", Expr::mul(A, x), Domain::equalsTo(b));
	M->solve();

        M->acceptedSolutionStatus(AccSolutionStatus::Anything);

        double tm = M->getSolverDoubleInfo("optimizerTime");
        double opt = M->getSolverDoubleInfo("intpntPrimalObj");

        auto Ax = *(M->getConstraint("linEq")->level());
        double Axbnrm = 0, bnrm = 0;

        for (int i = 0; i < bref.size(); i++) {
            Axbnrm += (Ax[i] - (*b)[i]) * (Ax[i] - (*b)[i]);
            bnrm += (*b)[i] * (*b)[i];
        }

        auto z = *(x->dual());
        auto y = *(M->getConstraint("linEq")->dual());
        std::vector<double> Aty(MatrixA(problemData).cols(), 0);

        for (int i = 0; i < MatrixA(problemData).nonZeros(); i++) {
            Aty[(*subj)[i]] -= y[(*subi)[i]] * ((*vals)[i]);
        }
        ind = 1;
        for (auto &psdMatrixSize : psdMatrixSizes(problemData)) {
            int cumsum = psdMatrixSize;
            int cnt = 0;
            for (int i = 1; i <= psdMatrixSize * (psdMatrixSize + 1) / 2; i++) {
                if (cnt == 0) { cnt = cumsum - 1; cumsum--; continue; }
                Aty[ind + i - 1] /= 2; cnt--;
            }
            ind += psdMatrixSize * (psdMatrixSize + 1) / 2;
        }

        double Atynrm = 0, cnrm = 1;
        for (int i = 0; i < MatrixA(problemData).cols(); i++) {
            if (i == 0) {
                Atynrm += (z[i] - Aty[i] - 1) * (z[i] - Aty[i] - 1);
            }
            else {
                Atynrm += (z[i] - Aty[i]) * (z[i] - Aty[i]);
            }
        }

        double pinf = sqrt(Axbnrm) / (1 + sqrt(bnrm));
        double dinf = sqrt(Atynrm) / 2;
        double gap = std::abs(M->getSolverDoubleInfo("intpntPrimalObj") - M->getSolverDoubleInfo("intpntDualObj")) / (1 + std::abs(M->getSolverDoubleInfo("intpntPrimalObj")) + std::abs(M->getSolverDoubleInfo("intpntDualObj")));
        double err = std::max({pinf, dinf, gap});
        int ite = M->getSolverIntInfo("intpntIter");
        if (verbose(problemData)) {
            std::cout << "time : " <<  tm << ", opt : " << opt << ", err :" << err  << ", ite : " << ite << std::endl;
        }
        //std::cout << "Solution : " << (*(x->level()))[0] << std::endl;

        std::vector<double> _y(y.size());
        for (int i = 0; i < y.size(); i++) _y[i] = y[i];
        return std::make_pair(opt, _y);
    }
}
