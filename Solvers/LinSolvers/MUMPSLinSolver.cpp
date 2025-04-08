#include "MUMPSLinSolver.h"

MUMPSLinSolver::MUMPSLinSolver(){
    MyName = "MUMPS";

    PetscCallThrow(KSPCreate(MPI_COMM_WORLD, &ksp));
    PetscCallThrow(KSPSetType(ksp, KSPPREONLY));
    PetscCallThrow(KSPSetErrorIfNotConverged(ksp, PETSC_TRUE));
    PetscCallThrow(KSPGetPC(ksp, &pc));
    PetscCallThrow(PCSetType(pc, PCLU));
    PetscCallThrow(PCFactorSetMatSolverType(pc, MATSOLVERMUMPS));
}

MUMPSLinSolver::~MUMPSLinSolver(){
    KSPDestroy(&ksp);
}

void MUMPSLinSolver::init(inputData& inputs, Mat& K){
    inputs.GetRequired(IOptionIdx, {"nonlinSolver","LinSolver","IOptionIdx"});
    inputs.GetRequired(IOptionVal, {"nonlinSolver","LinSolver","IOptionVals"});
}

void MUMPSLinSolver::solve(Mat& K, Vec& f, Vec& dx){
    PetscCallThrow(KSPSetOperators(ksp, K, K));
    PetscCallThrow(PCFactorSetUpMatSolverType(pc));

    Mat MumpsMat;
    PetscCallThrow(PCFactorGetMatrix(pc,&MumpsMat));
    for (size_t i = 0; i < IOptionIdx.size(); i++){
        PetscCallThrow(MatMumpsSetIcntl(MumpsMat, IOptionIdx[i], IOptionVal[i]));
    }  
    
    PetscCallThrow(KSPSetUp(ksp));
    PetscCallThrow(KSPSolve(ksp, f, dx));
}