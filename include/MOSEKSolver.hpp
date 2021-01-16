#ifndef __MOSEK_SOLVER_HPP__
#define __MOSEK_SOLVER_HPP__

#include "Solver.hpp"
#include "mosek/fusion.h"

namespace SPOPT {

    class MOSEKSolver : public Solver {
        public:
            // constructors
            MOSEKSolver(){};
            // Reads solver settings from the file `filename` and store it on `solverConfig`
            // when instanciating the Solver.
            // Internally, it calls a function LoadConfig.
            MOSEKSolver(std::string fileName);

            void LoadConfig(std::string fileName);
            std::pair<double, std::vector<double>> Solve(const ProblemData &problemData);
    };
}
#endif //__HSDE_SOLVER_HPP__