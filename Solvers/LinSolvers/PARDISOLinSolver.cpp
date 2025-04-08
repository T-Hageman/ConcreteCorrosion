#include "PARDISOLinSolver.h"

PARDISOLinSolver::PARDISOLinSolver(){
    MyName = "Pardiso";
    firstSolve = true;
}

PARDISOLinSolver::~PARDISOLinSolver(){

}

void PARDISOLinSolver::init(inputData& inputs, Mat& K){
    inputs.GetRequired(IOptionIdx, {"nonlinSolver","LinSolver","IOptionIdx"});
    inputs.GetRequired(IOptionVal, {"nonlinSolver","LinSolver","IOptionVals"});
}

void PARDISOLinSolver::solve(Mat& K, Vec& f, Vec& dx){
    if (firstSolve){
        PetscCallThrow(KSPCreate(MPI_COMM_WORLD, &ksp));
        PetscCallThrow(KSPSetType(ksp, KSPPREONLY));
        PetscCallThrow(KSPSetErrorIfNotConverged(ksp, PETSC_TRUE));
        PetscCallThrow(KSPGetPC(ksp, &pc));
        PetscCallThrow(PCSetType(pc, PCLU));
        PetscCallThrow(PCFactorSetMatSolverType(pc, MATSOLVERMKL_CPARDISO));
    }
    PetscCallThrow(KSPSetOperators(ksp, K, K));

    PetscCallThrow(PCFactorSetUpMatSolverType(pc));
    Mat PardisoMat;
    PetscCallThrow(PCFactorGetMatrix(pc,&PardisoMat));
    for (size_t i = 0; i < IOptionIdx.size(); i++){
        PetscCallThrow(MatMkl_PardisoSetCntl(PardisoMat, IOptionIdx[i], IOptionVal[i]));
    }  
    
    PetscCallThrow(KSPSetUp(ksp));
    PetscCallThrow(KSPSolve(ksp, f, dx));
    //KSPDestroy(&ksp);

    firstSolve = false;
}