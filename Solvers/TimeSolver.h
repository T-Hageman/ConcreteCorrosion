#ifndef TIMELINSOLVER_H
#define TIMELINSOLVER_H

#include "NonLinSolver.h"
#include "../InputsOutputs/OutputSaver.h"

class Restarter;

/// @brief Implements standard time-stepping scheme
class TimeSolver {
    public:
        TimeSolver(inputData& inputs, Physics& physics_in, NonLinSolver& NLSolver_in, Restarter& Restart_in);
        TimeSolver(inputData& inputs, Physics& physics_in, NonLinSolver& NLSolver_in, Restarter& Restart_in ,SaveDataFile& data);
        void Setup(inputData& inputs, Physics& physics_in, NonLinSolver& NLSolver_in, Restarter& Restart_in);
        void Save(SaveDataFile& data);
        ~TimeSolver();

        void Solve();
    private:
        PetscInt size, rank;
        size_t step;            //current timestep
        size_t outputN;         //output every N timesteps
        double t, dt, dt0, tmax;     //current, increment, initialization increment, and maximum time
        double dtGrow, dtMax;
        double dtScale;

        bool UseDtSteps, DynamicDT, RefineByData;
        std::string TimeRefineDataName; 
        double TimeRefineDataLim;
        std::vector<double> dtSteps, tSteps;

        Physics* physics;       //pointer to physics object
        NonLinSolver* NRSolver; //pointer to nonlinear solver object
        OutputDataSaver* Saver; //pointer to output-saving object
        Restarter* RestartSaver;//pointer to restart-saving object

        void Setdt();
};
#endif