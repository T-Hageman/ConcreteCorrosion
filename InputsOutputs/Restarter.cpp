#include "Restarter.h"
#include <filesystem>

/// @brief Initialize based on input files
/// @param inputs 
Restarter::Restarter(inputData& inputs){
    // Load input properties
    inputs.GetRequired(Restartable, {"Restart","Restartable"});
    inputs.GetRequired(SaveFreq, {"Restart","SaveFrequency"});
    inputs.GetRequired(nBackupsMax, {"Restart","BackUps"});
    nBackupsExist = 0;
    inputs.GetRequired(File_Prefix, {"Restart","Saveloc"});

    //check if backups already exist
    if (std::filesystem::exists(File_Prefix+".hdf5")){
        nBackupsExist = 1;
        for (size_t i = 0; i < nBackupsMax; i++){
            if (std::filesystem::exists(File_Prefix+"Old"+std::to_string(i+1)+".hdf5")){
                nBackupsExist+= 1;
            }
        }
    }

    // Determine if the current simulations can be started from restart files
    if (Restartable && nBackupsExist>0){
        Restartable = true;
    } else {
        Restartable = false;
    }
}

Restarter::~Restarter(){


}

/// @brief Loads simulation data to restart from files
/// @param inputs Input data object
void Restarter::Restart(inputData& inputs){
    SaveDataFile dataFile(File_Prefix+".hdf5", "read");

    dataFile.SetPrefix("Mesh");
    MyMesh = new Mesh(inputs, dataFile);

    dataFile.SetPrefix("Physics");
    MyPhysics = new Physics(inputs, *MyMesh, dataFile);

    dataFile.SetPrefix("NRSolver");
    NRSolver = new NonLinSolver(inputs, *MyPhysics, dataFile);

    dataFile.SetPrefix("TimeSolver");
    TimeStepper = new TimeSolver(inputs, *MyPhysics, *NRSolver, *this, dataFile);
}

/// @brief Register simulation componenets such that this class can save outputs for them
/// @param mesh_in //pointer to mesh object
/// @param physics_in //pointer to physics object
/// @param nonlinsolver_in //pointer to nonlinear solver
/// @param timesolver_in //pointer to time solver
void Restarter::Register(Mesh* mesh_in, Physics* physics_in, NonLinSolver* nonlinsolver_in, TimeSolver* timesolver_in){
    MyMesh = mesh_in;
    MyPhysics = physics_in;
    NRSolver = nonlinsolver_in;
    TimeStepper = timesolver_in;
}

/// @brief Saves a restart file
/// @param step Current time step
void Restarter::Save(size_t step){

    if (step%SaveFreq == 0){
        PetscMPIInt    size; PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&size));
        PetscMPIInt    rank; PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));

        //management of previously saved restart files
        if (rank == 0){
            if (std::filesystem::exists(File_Prefix+std::to_string(nBackupsMax)+".hdf5")){
                std::filesystem::remove(File_Prefix+std::to_string(nBackupsMax)+".hdf5");
            }
            for (int i = nBackupsMax-1; i > 0; i--){
                if (std::filesystem::exists(File_Prefix+std::to_string(i)+".hdf5")){
                    std::filesystem::rename(File_Prefix+std::to_string(i)+".hdf5", File_Prefix+std::to_string(i+1)+".hdf5");
                }
            }
            if (std::filesystem::exists(File_Prefix+".hdf5")){
                std::filesystem::rename(File_Prefix+".hdf5", File_Prefix+std::to_string(1)+".hdf5");
            }
        }
        PetscBarrier(NULL);

        // save to new file
        SaveDataFile dataFile(File_Prefix+".hdf5", "write");
        dataFile.SetPrefix("Mesh");
        MyMesh->Save(dataFile);

        dataFile.SetPrefix("Physics");
        MyPhysics->Save(dataFile);

        dataFile.SetPrefix("Physics");
        NRSolver->Save(dataFile);

        dataFile.SetPrefix("TimeSolver");
        TimeStepper->Save(dataFile);
    }
}

Mesh* Restarter::GetMesh(){
    return MyMesh;
}
Physics* Restarter::GetPhysics(){
    return MyPhysics;
}
NonLinSolver* Restarter::GetNonLinSolver(){
    return NRSolver;
}
TimeSolver* Restarter::GetTimeSolver(){
    return TimeStepper;
}