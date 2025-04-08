
static char help[] = "Corrosion simulation code.\n\n";
  
#include <chrono>

#include "InputsOutputs/inputData.h"
#include "mesh/mesh.h"
#include "utility/utility.h"
#include "Physics/physics.h"
#include "Solvers/NonLinSolver.h"
#include "Solvers/TimeSolver.h"
#include "InputsOutputs/Restarter.h"
 
int Simulation(char **args){
   char const* InputName = args[1];
   inputData inputFile(InputName);
   Logs.Init(&inputFile);
   Restarter Restart(inputFile);

   Mesh* MyMesh;
   Physics* MyPhysics;
   NonLinSolver* NRSolver;
   TimeSolver* TimeStepper;
   if (Restart.Restartable){
      Restart.Restart(inputFile);
      MyMesh = Restart.GetMesh();
      MyPhysics = Restart.GetPhysics();
      NRSolver = Restart.GetNonLinSolver();
      TimeStepper = Restart.GetTimeSolver();
   } else {
      MyMesh = new Mesh(inputFile);
      MyPhysics = new Physics(inputFile, *MyMesh);
      NRSolver = new NonLinSolver(inputFile, *MyPhysics);
      TimeStepper = new TimeSolver(inputFile, *MyPhysics, *NRSolver, Restart);
      Restart.Register(MyMesh, MyPhysics, NRSolver, TimeStepper);
   }

   TimeStepper->Solve();

   delete TimeStepper;
   delete NRSolver;
   delete MyPhysics;
   delete MyMesh;
   return 0; 
}  
   
int main(int argc,char **args)
{
   PetscMPIInt    size;
   PetscMPIInt    rank;
   std::chrono::high_resolution_clock::time_point start, end;

   start = std::chrono::high_resolution_clock::now();
   PetscInitialize(&argc,&args,(char*)0,help);
   PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&size));
   PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
   PetscSetFPTrap(PETSC_FP_TRAP_ON);

   Logs.PrintSingle("Starting Code\n",0);
   Logs.PrintEvery("Core "+std::to_string(rank)+"/"+std::to_string(size)+"\n",0);

   RegisterModels();
   Simulation(args);

   PetscCallThrow(PetscBarrier(NULL));
   PetscFinalize();

   end = std::chrono::high_resolution_clock::now();
   if (rank==0){
      int hours = std::chrono::duration_cast<std::chrono::hours>(end-start).count();
      int minutes = std::chrono::duration_cast<std::chrono::minutes>(end-start).count()-60*hours;
      int seconds = std::chrono::duration_cast<std::chrono::seconds>(end-start).count()-3600*hours-60*minutes;
      std::cout << "Runtime: " << hours << " Hours, " << minutes << " Minutes, " << seconds << " Seconds\n";
   }

   return 0;
}

