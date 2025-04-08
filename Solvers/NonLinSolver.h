#ifndef NONLINSOLVER_H
#define NONLINSOLVER_H

#include "../Physics/physics.h"
#include "../InputsOutputs/OutputSaver.h"
#include "LinSolvers/LinSolveRegister.h"

/// @brief Nonlinear, Newton-Raphson type, solver with ability to resolve equations in an iteratively staggered manner
class NonLinSolver {
    public:
        NonLinSolver(inputData& inputs, Physics & physics_in);
        NonLinSolver(inputData& inputs, Physics& physics_in, SaveDataFile& data);
        void Setup(inputData& inputs, Physics& physics_in);
        ~NonLinSolver();

        void Save(SaveDataFile& data);
        bool SolveStep();

        Physics* physics;   //pointer to underlying physics object
        std::vector<Vec> dx;    //state increment
        std::vector<BaseLinSolver*> LinSolvers;   //Linear solver to use

    private:
        bool CheckConvergence(size_t step);
        void LinLineSearch(size_t step);
        void HOLineSearch(size_t step);

        bool SolveStep_Standard(size_t step);

        PetscInt size, rank;

        size_t maxIt, max_outer_it; //maximum number of iterations and staggered loops to perform
        std::vector<double> tiny;    //absolute value below which the solution is considered converged, independent on normalisation
        size_t nSteps;  //numer of staggered steps within solving procedure
        std::vector<double> convCrit;   //state and force convergence criteria

        bool lineSearch;    //flag whether to perform a linear line-search
        std::vector<double> lineSearchLimits; //boundaries of line-search algorithm
        std::vector<std::vector<double>> Err0, Err, ErrIt0; //error after solving step 0 and current error, for each staggered step

        OutputDataSaver* Saver; //object used to save outputData (FOR DEBUGGING PURPOSES ONLY)

        std::vector<bool> BFGS;  //flag whether to use BFGS for each step
        void BFGSSetup(inputData& inputs, size_t step);
        bool SolveStep_BFGS(size_t step);
        std::vector<std::vector<Vec>> dX, dr;
        std::vector<Vec> q;
        std::vector<std::vector<double>> bfgs_rho, bfgs_a, bfgs_b;
        size_t BFGSReset;
};
#endif