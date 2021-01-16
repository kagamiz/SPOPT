#include "SCSSolver.hpp"

namespace SPOPT {

    SCSSolver::SCSSolver(std::string fileName)
    {
        LoadConfig(fileName);
    }

    void SCSSolver::LoadConfig(std::string fileName)
    {
        Solver::LoadConfig(fileName);
        SCSParam.alpha    = solverConfig["alpha"].as<double>(1.0);
    }

    void SCSSolver::SetUpFrom(const ProblemData &problemData)
    {
        d = (ScsData *)scs_calloc(1, sizeof(ScsData));
        d->m = MatrixA(problemData).cols();
        d->n = MatrixA(problemData).rows();
        d->b = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
        d->c = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
        d->b[0] = -1;
        for (int i = 1; i < d->m; i++) {
            d->b[i] = 0;
        }
        for (int i = 0; i < d->n; i++) {
            d->c[i] = VectorB(problemData)[i] / (VectorD(problemData)[i] * dualScaler(problemData) * scalingFactor(problemData));
        }

        scs_int Anz;
        d->A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
        d->A->m = d->m;
        d->A->n = d->n;
        d->A->p = (scs_int *)scs_calloc(d->A->n + 1, sizeof(scs_int));
        Anz = MatrixA(problemData).nonZeros();
        d->A->x = (scs_float *)scs_calloc(Anz, sizeof(scs_float));
        d->A->i = (scs_int *)scs_calloc(Anz, sizeof(scs_int));

        std::vector<int> AColPerm(MatrixA(problemData).cols());
        AColPerm[0] = 0;
        int ctr = 1;
        int IndOffset = 1;
        int freeConeSize = 1;
        for (int i = 0; i < psdMatrixSizes(problemData).size(); i++) {
            IndOffset += psdMatrixSizes(problemData)[i] * (psdMatrixSizes(problemData)[i] + 1) / 2;
        }
        for (int i = 0; i < symmetricMatrixSizes(problemData).size(); i++) {
            int symSize = symmetricMatrixSizes(problemData)[i] * (symmetricMatrixSizes(problemData)[i] + 1) / 2;
            for (int j = 0; j < symSize; j++) {
                AColPerm[ctr] = IndOffset;
                IndOffset++; ctr++;
            }
            freeConeSize += symSize;
        }
        IndOffset = 1;
        for (int i = 0; i < psdMatrixSizes(problemData).size(); i++) {
            int psdSize = psdMatrixSizes(problemData)[i] * (psdMatrixSizes(problemData)[i] + 1) / 2;
            for (int j = 0; j < psdSize; j++) {
                AColPerm[ctr] = IndOffset;
                IndOffset++; ctr++;
            }
        }
        Eigen::SparseMatrix<double> A2 = MatrixA(problemData);
        int ind2 = 0;
        for (int i = 0; i < A2.outerSize(); i++) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A2, i); it; ++it) {
                while (ind2 < sq2Cols(problemData).size() && sq2Cols(problemData)[ind2] < it.col()) ind2++;
                double scalingCoef = (ind2 < sq2Cols(problemData).size() && sq2Cols(problemData)[ind2] == it.col() && it.col() >= IndOffset) ? sqrt(0.5) : 1;
                it.valueRef() = -it.value() / (scalingCoef * VectorD(problemData)[it.row()] * VectorE(problemData)[it.col()] * scalingFactor(problemData));
            }
        }

        std::vector<Eigen::Triplet<double>> tripletList;
        for (int i = 0; i < A2.outerSize(); i++) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(A2, AColPerm[i]); it; ++it) {
                tripletList.emplace_back(i, it.row(), it.value());
            }
        }
        Eigen::SparseMatrix<double> At;
        At.resize(MatrixA(problemData).cols(), MatrixA(problemData).rows());
        At.setFromTriplets(tripletList.begin(), tripletList.end());
        
        int ind = 0;
        for (int i = 0; i < d->A->n; i++) {
            d->A->p[i] = ind;
            for (Eigen::SparseMatrix<double>::InnerIterator it(At, i); it; ++it) {
                d->A->i[ind] = it.row();
                d->A->x[ind] = it.value();
                ind++;
            }
        }
        d->A->p[d->A->n] = ind;

        d->stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
        d->stgs->normalize             = enableScaling(problemData);
        d->stgs->scale                 = (enableScaling(problemData) ? 5 : 1);
        d->stgs->rho_x                 = RHO_X;
        d->stgs->max_iters             = basicParam.iterationLimit;
        d->stgs->eps                   = std::min({basicParam.primalTolerance, basicParam.dualTolerance, basicParam.gapTolerance});
        d->stgs->alpha                 = SCSParam.alpha;
        d->stgs->cg_rate               = CG_RATE;
        d->stgs->verbose               = verbose(problemData);
        d->stgs->warm_start            = false;
        d->stgs->acceleration_lookback = (basicParam.enableAndersonAcceleration ? basicParam.AAMemoryLength : 0);

        k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
        k->f = freeConeSize;
        k->l = 0;
        k->qsize = 0;
        k->q = (scs_int *)scs_calloc(k->qsize, sizeof(scs_int));
        k->ssize = psdMatrixSizes(problemData).size();
        k->s = (scs_int *)scs_calloc(k->ssize, sizeof(scs_int));
        for (int i = 0; i < k->ssize; i++) {
            k->s[i] = psdMatrixSizes(problemData)[i];
        }
        k->ep = 0;
        k->ed = 0;
        k->psize = 0;
        k->p = (scs_float *)scs_calloc(k->psize, sizeof(scs_float));
    }

    std::pair<double, std::vector<double>> SCSSolver::Solve(const ProblemData &problemData)
    {
        SetUpFrom(problemData);
        sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
        scs(d, k, sol, &info);

        if (verbose(problemData)) {
            std::cout << std::fixed << std::setprecision(2) << "time : " <<  (info.setup_time + info.solve_time) / 1000.0 << std::scientific << std::setprecision(2) << ", opt : " << info.pobj << ", err :" << std::max({info.res_pri, info.res_dual, info.rel_gap})  << ", ite : " << info.iter << std::endl;
            std::cout << std::resetiosflags(std::ios_base::floatfield); 
        }
        std::vector<double> y(MatrixA(problemData).rows());
        for (int i = 0; i < y.size(); i++) {
            y[i] = sol->x[i];
        }

        scs_free_data(d, k);
        scs_free_sol(sol);

        return std::make_pair(info.pobj, y);
    }   
}
