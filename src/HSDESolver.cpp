#include "HSDESolver.hpp"

namespace SPOPT {

    HSDESolver::HSDESolver(std::string fileName)
    {
        LoadConfig(fileName);
    }

    void HSDESolver::LoadConfig(std::string fileName)
    {
        Solver::LoadConfig(fileName);
        hsdeParam.alpha    = solverConfig["alpha"].as<double>(1.0);
    }

    void HSDESolver::SetUpFrom(const ProblemData &problemData)
    {
        At = MatrixA(problemData).transpose();
        Eigen::SparseMatrix<double> IAAtFull = MatrixA(problemData) * At;

        for (int i = 0; i < MatrixA(problemData).rows(); i++) {
            IAAtFull.coeffRef(i, i) += 1;
        }

        if (problemData.IsUnconstrained()) {
            sparseIAAt = IAAtFull.diagonal();
        }
        else {
            IAAt.compute(IAAtFull);
        }

        zeta = Eigen::VectorXd::Zero(MatrixA(problemData).cols() + MatrixA(problemData).rows());
        zeta.head(MatrixA(problemData).cols()) = -VectorC(problemData);
        zeta.tail(MatrixA(problemData).rows()) = VectorB(problemData);
        MInvZeta = CalcMInv(problemData, zeta);
        zetaDotMInvZeta = zeta.dot(MInvZeta);
    }

    Eigen::VectorXd HSDESolver::ConstructInitialPoint(const ProblemData &problemData)
    {
        int vecLength = MatrixA(problemData).cols() + MatrixA(problemData).rows() + 1;
        Eigen::VectorXd u = Eigen::VectorXd::Zero(vecLength);
        Eigen::VectorXd v = Eigen::VectorXd::Zero(vecLength);

        u(vecLength - 1) = 1;
        v(vecLength - 1) = 1;

        Eigen::VectorXd ret(2 * vecLength);
        ret << u, v;
        return ret;
    }

    Eigen::VectorXd HSDESolver::CalcMInv(const ProblemData &problemData, Eigen::VectorXd &v)
    {
        int rowNumOfA = MatrixA(problemData).rows();
        int colNumOfA = MatrixA(problemData).cols();
        
        Eigen::VectorXd ret = Eigen::VectorXd::Zero(colNumOfA + rowNumOfA);
        ret.tail(rowNumOfA) = MatrixA(problemData) * v.head(colNumOfA) + v.tail(rowNumOfA);
        if (problemData.IsUnconstrained()) {
            for (int i = 0; i < rowNumOfA; i++) {
                ret(colNumOfA + i) /= sparseIAAt(i);
            }
        }
        else {
            ret.tail(rowNumOfA) = IAAt.solve(ret.tail(rowNumOfA));
        }
        ret.head(colNumOfA) = v.head(colNumOfA) - At * ret.tail(rowNumOfA);
        return ret;
    }

    Eigen::VectorXd HSDESolver::ApplyFixedPointFunction(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        int vecLength = MatrixA(problemData).cols() + MatrixA(problemData).rows() + 1;
        Eigen::VectorXd curU = v.head(vecLength);
        Eigen::VectorXd curV = v.tail(vecLength);
        Eigen::VectorXd tmpVector(vecLength - 1);
        Eigen::VectorXd uHat(vecLength);

        double coef = curU(vecLength - 1) + curV(vecLength - 1);
        tmpVector = (curU + curV).head(vecLength - 1) - coef * zeta;
        tmpVector = CalcMInv(problemData, tmpVector);
        tmpVector = tmpVector - (zeta.dot(tmpVector) / (1 + zetaDotMInvZeta)) * MInvZeta;
        uHat.head(vecLength - 1) = tmpVector;
        uHat(vecLength - 1) = coef - VectorC(problemData).dot(tmpVector.head(MatrixA(problemData).cols())) + VectorB(problemData).dot(tmpVector.tail(MatrixA(problemData).rows()));

        uHat = hsdeParam.alpha * uHat + (1 - hsdeParam.alpha) * curU;

        Eigen::VectorXd newU = uHat - curV;
        if (enableLowerBound(problemData)) {
            newU(0) = std::max(newU(0), 0.0);
        }
        if (enableUpperBound(problemData)) {
            newU(enableLowerBound(problemData)) = std::max(newU(enableLowerBound(problemData)), 0.0);
        }
        int leftmostPosition = 1 + enableLowerBound(problemData) + enableUpperBound(problemData);
        for (auto psdMatrixSize : psdMatrixSizes(problemData)) {
            Eigen::MatrixXd tmpMatrix(psdMatrixSize, psdMatrixSize);

            for (int i = 0; i < psdMatrixSize; i++) {
                for (int j = 0; j < psdMatrixSize; j++) {
                    tmpMatrix(i, j) = newU(leftmostPosition + i * psdMatrixSize + j);
                    if (i > j) {
                        double p = tmpMatrix(i, j), q = tmpMatrix(j, i);
                        tmpMatrix(i, j) = tmpMatrix(j, i) = (p + q) / 2;
                    }
                }
            }

            Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(tmpMatrix);
            if (eigenSolver.info() != Eigen::Success) {
                std::cerr << "[ERROR] Eigen decomposition failed!" << std::endl;
                exit(1);
            }

            Eigen::MatrixXd tmpMatrix2 = Eigen::MatrixXd::Zero(psdMatrixSize, psdMatrixSize);
            for (int i = psdMatrixSize - 1; i >= 0; i--) {
                double lambda = eigenSolver.eigenvalues()(i);
                Eigen::VectorXd ivec = eigenSolver.eigenvectors().col(i);
                if (lambda > 0) {
                    tmpMatrix2 += lambda * ivec * ivec.transpose();
                }
                else {
                    break;
                }
            }

            for (int i = 0; i < psdMatrixSize; i++) {
                for (int j = 0; j < psdMatrixSize; j++) {
                    newU(leftmostPosition + i * psdMatrixSize + j) = tmpMatrix2(i, j);
                }
            }
            leftmostPosition += psdMatrixSize * psdMatrixSize;
        }
        if (newU(vecLength - 1) < 0) newU(vecLength - 1) = 0;

        Eigen::VectorXd newV = curV - uHat + newU;

        Eigen::VectorXd ret(2 * vecLength);
        ret << newU, newV;
        return ret;
    }

    double HSDESolver::GetPrimalInfeasibility(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        int vecLength = MatrixA(problemData).rows() + MatrixA(problemData).cols() + 1;
        Eigen::VectorXd vecU = v.head(vecLength), vecV = v.tail(vecLength);

        if (vecU(vecLength - 1) <= 0) {
            return std::numeric_limits<double>::max();
        }

        Eigen::VectorXd x = vecU.head(MatrixA(problemData).cols()) / vecU(vecLength - 1);
        Eigen::VectorXd primalResidual = MatrixA(problemData) * x - VectorB(problemData);

        for (int i = 0; i < MatrixA(problemData).rows(); i++) {
            primalResidual(i) /= (dualScaler(problemData) * VectorD(problemData)(i));
        }

        return primalResidual.norm() / (1 + originalBNorm(problemData));
    }

    double HSDESolver::GetDualInfeasibility(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        int vecLength = MatrixA(problemData).rows() + MatrixA(problemData).cols() + 1;
        Eigen::VectorXd vecU = v.head(vecLength), vecV = v.tail(vecLength);

        if (vecU(vecLength - 1) <= 0) {
            return std::numeric_limits<double>::max();
        }

        Eigen::VectorXd y = vecU.segment(MatrixA(problemData).cols(), MatrixA(problemData).rows()) / vecU(vecLength - 1);
        Eigen::VectorXd z = vecV.head(MatrixA(problemData).cols()) / vecU(vecLength - 1);

        Eigen::VectorXd dualResidual = z - At * y + VectorC(problemData);

        for (int i = 0; i < MatrixA(problemData).cols(); i++) {
            dualResidual(i) /= (primalScaler(problemData) * VectorE(problemData)(i));
        }

        return dualResidual.norm() / (1 + originalCNorm(problemData));
    }

    double HSDESolver::GetPrimalObjValue(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        int vecLength = MatrixA(problemData).rows() + MatrixA(problemData).cols() + 1;
        Eigen::VectorXd vecU = v.head(vecLength);

        if (vecU(vecLength - 1) <= 0) {
            return std::numeric_limits<double>::max();
        }

        Eigen::VectorXd x = vecU.head(MatrixA(problemData).cols()) / vecU(vecLength - 1);

        return VectorC(problemData).dot(x) / (primalScaler(problemData) * dualScaler(problemData));
    }

    double HSDESolver::GetDualObjValue(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        int vecLength = MatrixA(problemData).rows() + MatrixA(problemData).cols() + 1;
        Eigen::VectorXd vecU = v.head(vecLength);

        if (vecU(vecLength - 1) <= 0) {
            return std::numeric_limits<double>::max();
        }

        Eigen::VectorXd y = vecU.segment(MatrixA(problemData).cols(), MatrixA(problemData).rows()) / vecU(vecLength - 1);

        return VectorB(problemData).dot(y) / (primalScaler(problemData) * dualScaler(problemData));
    }

    double HSDESolver::GetGap(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        int vecLength = MatrixA(problemData).rows() + MatrixA(problemData).cols() + 1;
        Eigen::VectorXd vecU = v.head(vecLength), vecV = v.tail(vecLength);

        if (vecU(vecLength - 1) <= 0) {
            return std::numeric_limits<double>::max();
        }

        double primalObjVal = GetPrimalObjValue(problemData, v);
        double dualObjVal   = GetDualObjValue(problemData, v);

        double slacknessInfeasibility = dualObjVal - primalObjVal;

        return abs(slacknessInfeasibility) / (1 + abs(primalObjVal) + abs(dualObjVal));
    }

    void HSDESolver::UpdateParameter(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
    }
}