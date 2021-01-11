#ifndef __PROBLEM_DATA_HPP__
#define __PROBLEM_DATA_HPP__

#include <yaml-cpp/yaml.h>

#include "LinearAlgebra.hpp"
#include "Polynomial.hpp"
#ifdef __BUILD_WITH_MATLAB__
#include "MATLAB.hpp"
#endif //__BUILD_WITH_MATLAB__

namespace SPOPT {

    class Solver;
    class Polynomial;

    // stores the information of the POP problem :
    //
    // min {x \in K} f(x)
    // where K = {x | g_i(x) >= 0 (i = 1, ..., numIneq), h_j(x) = 0 (j = 1, ..., numEq)}
    //
    // also, when the relaxation degree "d" of the Lasserre hierarchy is given in configuration file,
    // this class will store the information of corresponding SDP's data matrix and vectors.

    class ProblemData {
        public:
            // constructors
            ProblemData(){};
            // Reads solver settings from the file `filename` and store the information on private members.
            // when instanciating the ProblemData.
            // Internally, it calls a function LoadConfig.
            ProblemData(std::string fileName);

            // Reads solver settings from the file `filename` and store the information on private members.
            // TODO: error handling
            void LoadConfig(std::string fileName);

            void ConstructSDP();
            void ShowPOP();

            // check the function value and the gradient at given point
            void Analyze(std::vector<double> &x);

            // Output for other solvers
            void OutputJuliaFile(std::string fileName = "");
            // Output the Matrix A, b, c, and cone information
            void OutputMatFile(std::string fileName = "outfile.mat");

            // returns whether the problem is unconstrained problem (i.e. K = \mathbb{R}^n) or not.
            bool IsUnconstrained() const;

            // returns approximate point that satisfies moment conditions by
            // an algorithm proposed by Henrion & Lasserre (2006).
            std::vector<std::vector<double>> ExtractSolutionsFrom(std::vector<double> &v);
            
            friend Solver;

        private:
            /* Declaration of EnumType */
            enum ConstraintType {
                EqualityConstraint,
                InequalityConstraint,
                BoundConstraint
            };
            /* Data to be loaded from the function `LoadConfig`. */

            // Polynomials that characterize the POP
            Polynomial objectiveFunction;
            std::vector<Polynomial> originalInequalityConstraints;
            std::vector<Polynomial> originalEqualityConstraints;

            std::vector<std::vector<int>> originalIndexSets;
            std::vector<std::vector<int>> originalJunctionTree;

            // Options
            bool enableScaling;
            bool enableGradientConstraint;
            bool enableGradientConstraintType2;
            bool enableLowerBound;
            double lowerBoundConstant;
            bool enableUpperBound;
            double upperBoundConstant;
            double scalingFactor;

            /* Other data to be constructed from the function `ConstructSDP` */

            int hierarchyDegree;

            std::vector<Polynomial> convertedInequalityConstraints;
            std::vector<int>        groupIDOfConvertedInequalityConstraints;
            std::vector<Polynomial> convertedEqualityConstraints;
            std::vector<int>        groupIDOfConvertedEqualityConstraints;

            std::vector<std::vector<int>> convertedIndexSets;
            std::vector<std::vector<Term>> convertedTermSets;

            // Matrix and vector data used in SDP
            Eigen::SparseMatrix<double> A;
            Eigen::SparseVector<double> c;
            Eigen::VectorXd b;

            std::vector<int> sq2Cols;

            double originalBNorm, originalCNorm;

            // Scaling data used for boosting ADMM
            // calculated in function `ConstructSDP` when the option `enableScaling` is set to `true`.
            Eigen::VectorXd D, E;       // Scaling matrix
            double primalScaler;
            double dualScaler;

            std::vector<int> psdMatrixSizes, symmetricMatrixSizes; 

            /* Private Member Functions */

            void _AddGradientConstraints();
            void _ConstructNewConstraints(); // make new constraints by traversing originalJunctionTree
            void _TraverseTree(int v, int p, int &newVariableIndex,
                              std::vector<std::vector<Monomial>> &objectiveMonomials,
                              std::vector<std::vector<int>> &objectiveIDs,
                              std::vector<std::vector<std::vector<Monomial>>> &constraintMonomials,
                              std::vector<std::vector<int>> &constraintIDs,
                              std::vector<ConstraintType>   &constraintTypes,
                              std::vector<bool> &visited,
                              int &variableOrder);
            
            void _ConstructTermMap();
            void _GenerateTerms(std::vector<int> &inds, std::vector<Term> &terms, Term &tmp, int pos, int left);

            void _ConstructVectorB();
            void _ConstructMatrixA();
            void _ConstructVectorC();
            void _ConstructScalingData();
            
            // stores elimination ordering of the member of the `convertedJunctionTree`.
            std::map<int, int> variableOrderMap;

            // transforms term to an integer ordered by perfect elimination ordering of variables
            std::map<Term, int> termToInteger;

            // auxiliary function used to extract an approximate solutions from truncated moments
            void _GenerateSolutions(int cur, std::vector<std::vector<std::vector<double>>> &solutionCandidates, std::vector<bool> &done, std::vector<double> &tmp, std::vector<std::vector<double>> &answers);
};
}
#endif //__PROBLEM_DATA_HPP__