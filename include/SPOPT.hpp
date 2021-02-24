#ifndef __SPOPT_HPP__
#define __SPOPT_HPP__

#include "StandardLibraries.hpp"
#include "Polynomial.hpp"
#include "ProblemData.hpp"
#include "DualLagrangianSolver.hpp"
#include "HSDESolver.hpp"
#ifdef __BUILD_WITH_MOSEK__
#include "MOSEKSolver.hpp"
#endif
#ifdef __BUILD_WITH_SCS__
#include "SCSSolver.hpp"
#endif
#endif //__SPOPT_HPP__