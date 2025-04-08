#ifndef RESTARTER_H
#define RESTARTER_H

#include <petsc.h>
#include <fstream>
#include <iostream>

#include "SaveData.h"
#include "inputData.h"
#include "../Solvers/NonLinSolver.h"
#include "../Solvers/TimeSolver.h"

/// @brief Class to aid in saving restart files, and restarting from these files
class Restarter {
    public:
        Restarter(inputData& inputs);
        ~Restarter();

        void Restart(inputData& inputs);
        void Register(Mesh* mesh_in, Physics* physics_in, NonLinSolver* nonlinsolver_in, TimeSolver* timesolver_in);
        void Save(size_t step);

        Mesh* GetMesh();
        Physics* GetPhysics();
        NonLinSolver* GetNonLinSolver();
        TimeSolver* GetTimeSolver();

        bool Restartable;   //indicates whether the current simulation has pre-existing restart files from which to continue

    private:
        size_t SaveFreq;    //Frequency to save output files (every n timesteps)
        size_t nBackupsMax, nBackupsExist; //Number of backup files to keep before overwriting oldest

        std::string File_Prefix;    //prefix/path for restart files

        Mesh* MyMesh;
        Physics* MyPhysics;
        NonLinSolver* NRSolver;
        TimeSolver* TimeStepper;
};

#endif