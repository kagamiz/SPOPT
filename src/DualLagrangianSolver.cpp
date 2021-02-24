#include "DualLagrangianSolver.hpp"

namespace SPOPT {
    static const double sq2 = sqrt(2);

    DualLagrangianSolver::DualLagrangianSolver(std::string fileName)
    {
        LoadConfig(fileName);
    }

    void DualLagrangianSolver::LoadConfig(std::string fileName)
    {
        Solver::LoadConfig(fileName);
        dualParam.mu    = solverConfig["mu"].as<double>(0.1);
        dualParam.rho   = solverConfig["rho"].as<double>(1.6);
        dualParam.gamma = solverConfig["gamma"].as<double>(0.5);
        dualParam.h1    = solverConfig["h1"].as<int>(20);
        dualParam.h2    = solverConfig["h2"].as<int>(150);
        dualParam.h3    = solverConfig["h3"].as<int>(300);
        dualParam.h4    = solverConfig["h4"].as<int>(50);
    }

    void DualLagrangianSolver::SetUpFrom(const ProblemData &problemData)
    {
        At = MatrixA(problemData).transpose();
        Eigen::SparseMatrix<double> AAtFull = MatrixA(problemData) * At;

        if (problemData.IsUnconstrained()) {
            sparseAAt = AAtFull.diagonal();
        }
        else {
            AAt.compute(AAtFull);
        }

        it_pinf = it_dinf = 0;
    }

    Eigen::VectorXd DualLagrangianSolver::ConstructInitialPoint(const ProblemData &problemData)
    {
        Eigen::VectorXd x = Eigen::VectorXd::Zero(MatrixA(problemData).cols());
        Eigen::VectorXd y = Eigen::VectorXd::Zero(MatrixA(problemData).rows());
        Eigen::VectorXd z = Eigen::VectorXd::Zero(MatrixA(problemData).cols());

        int leftmostPosition = 1;
        for (auto psdMatrixSize : psdMatrixSizes(problemData)) {
            int acc = 0;
            for (int i = 0; i < psdMatrixSize; i++) {
                x(leftmostPosition + acc) = 1;
                acc += psdMatrixSize - i;
            }
            leftmostPosition += psdMatrixSize * (psdMatrixSize + 1) / 2;
        }

        Eigen::VectorXd v(x.size() + y.size() + z.size());
        v << x, y, z;
        return v;
    }

    Eigen::VectorXd DualLagrangianSolver::CalcInvAAt(const ProblemData &problemData, Eigen::VectorXd &v)
    {
        Eigen::VectorXd retV(v.size());

        if (problemData.IsUnconstrained()) {
            for (int i = 0; i < v.size(); i++) {
                retV(i) = v(i) / sparseAAt(i);
            }
        }
        else {
            retV = AAt.solve(v);
        }

        return retV;
    }

    Eigen::VectorXd DualLagrangianSolver::ApplyFixedPointFunction(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        Eigen::VectorXd curX = v.head(MatrixA(problemData).cols());
        Eigen::VectorXd curY = v.segment(MatrixA(problemData).cols(), MatrixA(problemData).rows());
        Eigen::VectorXd curZ = v.tail(MatrixA(problemData).cols());

        Eigen::VectorXd newY = dualParam.mu * (VectorB(problemData) - MatrixA(problemData) * curX) - MatrixA(problemData) * curZ;
        newY -= MatrixA(problemData) * VectorC(problemData);
        newY = -CalcInvAAt(problemData, newY);

        Eigen::VectorXd tmpZ = At * newY - dualParam.mu * curX - VectorC(problemData);
        Eigen::VectorXd newZ = Eigen::VectorXd::Zero(MatrixA(problemData).cols());

        int leftmostPosition = 1;
        for (auto psdMatrixSize : psdMatrixSizes(problemData)) {
            Eigen::MatrixXd tmpMatrix(psdMatrixSize, psdMatrixSize);

            int ptr = 0;
            for (int i = 0; i < psdMatrixSize; i++) {
                for (int j = i; j < psdMatrixSize; j++) {
                    tmpMatrix(i, j) = tmpZ(leftmostPosition + ptr) * (i != j ? sq2 / 2 : 1);
                    if (i != j) {
                        tmpMatrix(j, i) = tmpMatrix(i, j);
                    }
                    ptr++;
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

            ptr = 0;
            for (int i = 0; i < psdMatrixSize; i++) {
                for (int j = i; j < psdMatrixSize; j++) {
                    newZ(leftmostPosition + ptr) = tmpMatrix2(i, j) * (i != j ? sq2 : 1);
                    ptr++;
                }
            }
            leftmostPosition += psdMatrixSize * (psdMatrixSize + 1) / 2;
        }

        Eigen::VectorXd newX = curX + (newZ - At * newY + VectorC(problemData)) / dualParam.mu;
        Eigen::VectorXd newV(v.size());
        newV << newX, newY, newZ;
        return newV;
    }

    double DualLagrangianSolver::GetPrimalInfeasibility(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        Eigen::VectorXd x = v.head(MatrixA(problemData).cols());
        Eigen::VectorXd primalResidual = MatrixA(problemData) * x - VectorB(problemData);

        for (int i = 0; i < MatrixA(problemData).rows(); i++) {
            primalResidual(i) /= (scalingFactor(problemData) * dualScaler(problemData) * VectorD(problemData)(i));
        }

        return primalResidual.norm() / (1 + originalBNorm(problemData));
    }

    double DualLagrangianSolver::GetDualInfeasibility(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        Eigen::VectorXd y = v.segment(MatrixA(problemData).cols(), MatrixA(problemData).rows());
        Eigen::VectorXd z = v.tail(MatrixA(problemData).cols());

        Eigen::VectorXd dualResidual = z - At * y + VectorC(problemData);

        for (int i = 0; i < MatrixA(problemData).cols(); i++) {
            dualResidual(i) /= (scalingFactor(problemData) * primalScaler(problemData) * VectorE(problemData)(i));
        }

        return dualResidual.norm() / (1 + originalCNorm(problemData));
    }

    double DualLagrangianSolver::GetPrimalObjValue(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        Eigen::VectorXd x = v.head(MatrixA(problemData).cols());
        return VectorC(problemData).dot(x) / (scalingFactor(problemData) * primalScaler(problemData) * dualScaler(problemData));
    }

    double DualLagrangianSolver::GetDualObjValue(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        Eigen::VectorXd y = v.segment(MatrixA(problemData).cols(), MatrixA(problemData).rows());
        return VectorB(problemData).dot(y) / (scalingFactor(problemData) * primalScaler(problemData) * dualScaler(problemData));
    }

    double DualLagrangianSolver::GetGap(const ProblemData &problemData, const Eigen::VectorXd &v)
    {

        double primalObjVal = GetPrimalObjValue(problemData, v);
        double dualObjVal   = GetDualObjValue(problemData, v);

        double slacknessInfeasibility = dualObjVal - primalObjVal;

        return abs(slacknessInfeasibility) / (1 + abs(primalObjVal) + abs(dualObjVal));
    }

    void DualLagrangianSolver::UpdateParameter(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        double pinf = GetPrimalInfeasibility(problemData, v);
        double dinf = GetDualInfeasibility(problemData, v);

        if (pinf <= dinf) {
            it_pinf++;
            it_dinf = 0;
            if (it_pinf >= dualParam.h4) {
                dualParam.mu = std::clamp(dualParam.gamma * dualParam.mu, dualParam.muMin, dualParam.muMax);
                it_pinf = 0;
            }
        }
        else {
            it_dinf++;
            it_pinf = 0;
            if (it_dinf >= dualParam.h4) {
                dualParam.mu = std::clamp((1 / dualParam.gamma) * dualParam.mu, dualParam.muMin, dualParam.muMax);
                it_dinf = 0;
            }
        }
    }

    std::vector<double> DualLagrangianSolver::GetDualVariable(const ProblemData &problemData, const Eigen::VectorXd &v)
    {
        Eigen::VectorXd y = v.segment(MatrixA(problemData).cols(), MatrixA(problemData).rows());
        for (int i = 0; i < MatrixA(problemData).rows(); i++) {
            y(i) *= (dualScaler(problemData) * VectorD(problemData)[i]);
        }
        return std::vector<double>(y.data(), y.data() + y.size());
    }
}