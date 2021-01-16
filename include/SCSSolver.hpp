#ifndef __SCS_SOLVER_HPP__
#define __SCS_SOLVER_HPP__

#include "Solver.hpp"
#include "scs/scs.h"
#include "scs/util.h"
#include "linsys/amatrix.h"

namespace SPOPT {

    class SCSParameter {
        public:
            double alpha;
    };

    class SCSSolver : public Solver {
        public:
            // constructors
            SCSSolver(){};
            // Reads solver settings from the file `filename` and store it on `solverConfig`
            // when instanciating the Solver.
            // Internally, it calls a function LoadConfig.
            SCSSolver(std::string fileName);

            void LoadConfig(std::string fileName);
            std::pair<double, std::vector<double>> Solve(const ProblemData &problemData);

        private:
            void SetUpFrom(const ProblemData &problemData);

            ScsData *d;
            ScsCone *k;
            ScsInfo info;
            ScsSolution *sol;

            SCSParameter SCSParam;
    };
}
#endif //__HSDE_SOLVER_HPP__