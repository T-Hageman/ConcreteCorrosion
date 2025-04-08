#include "GMRESLinSolver.h"

#include "../../utility/utility.h"

GMRESLinSolver::GMRESLinSolver(){
    MyName = "GMRES";
    firstSolve = true;
}

GMRESLinSolver::~GMRESLinSolver(){

}

void GMRESLinSolver::init(inputData& inputs, Mat& K){
    inputs.GetRequired(IOptionIdx, {"nonlinSolver","LinSolver","IOptionIdx"});
    inputs.GetRequired(IOptionVal, {"nonlinSolver","LinSolver","IOptionVals"});
}

void GMRESLinSolver::solve(Mat& K, Vec& f, Vec& dx){
    if (firstSolve){
        PetscCallThrow(KSPCreate(MPI_COMM_WORLD, &ksp));
        PetscCallThrow(KSPSetType(ksp, KSPGMRES));
        PetscCallThrow(KSPSetErrorIfNotConverged(ksp, PETSC_FALSE));
        PetscCallThrow(KSPSetTolerances(ksp, 1.0e-9, 1.0e-15, PETSC_DEFAULT, PETSC_DEFAULT));
        PetscCallThrow(KSPGMRESSetRestart(ksp, 100));

        if (Logs.info_Level>3){
            PetscViewerAndFormat *vf;
            PetscCallThrow(PetscViewerAndFormatCreate(PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_DEFAULT, &vf));
            PetscCallThrow(KSPMonitorSet(ksp, (PetscErrorCode(*)(KSP, PetscInt, PetscReal, void *))KSPMonitorResidual, vf, (PetscErrorCode(*)(void **))PetscViewerAndFormatDestroy));
        }

        PetscCallThrow(KSPSetDiagonalScale(ksp, PETSC_TRUE));
        PetscCallThrow(KSPSetDiagonalScaleFix(ksp, PETSC_TRUE));

        if (true){
            PetscCallThrow(KSPGetPC(ksp, &pc));
            PetscCallThrow(PCSetType(pc, PCASM));
            PCASMSetOverlap(pc, 1);

            PetscCallThrow(KSPSetOperators(ksp, K, K));
            PetscCallThrow(KSPSetUp(ksp));

            KSP      *subksp;        /* array of KSP contexts for local subblocks */
            PetscInt  nlocal, first; /* number of local subblocks, first local subblock */
            PC        subpc;         /* PC context for subblock */

            PetscCallThrow(PCASMGetSubKSP(pc, &nlocal, &first, &subksp));
            for (int i = 0; i < nlocal; i++) {
                PetscCallThrow(KSPGetPC(subksp[i], &subpc));

                if (false){
                    PetscCallThrow(KSPSetType(subksp[i], KSPPREONLY));
                    PetscCallThrow(PCSetType(subpc, PCLU));
                    PetscCallThrow(PCFactorSetMatSolverType(subpc, MATSOLVERMKL_PARDISO));

                    Mat PardisoMat;
                    PetscCallThrow(PCFactorGetMatrix(subpc,&PardisoMat));
                    for (size_t i = 0; i < IOptionIdx.size(); i++){
                        PetscCallThrow(MatMkl_PardisoSetCntl(PardisoMat, IOptionIdx[i], IOptionVal[i]));
                    }  
                } else {
                    // PetscCallThrow(PCSetType(subpc, PCILU));
                    //PetscCallThrow(KSPSetType(subksp[i], KSPGMRES));
                    //PetscCallThrow(KSPSetTolerances(subksp[i], 1.e-9, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT));

                    PetscCallThrow(KSPSetType(subksp[i], KSPPREONLY));
                    PetscCallThrow(PCSetType(subpc, PCLU));
                    PetscCallThrow(PCFactorSetMatSolverType(subpc, MATSOLVERMUMPS));
                }
            }

            PetscCallThrow(PCSetReusePreconditioner(pc, PETSC_FALSE));
        

        } else {
            PetscCallThrow(KSPGetPC(ksp, &pc));
            PetscCallThrow(PCSetType(pc, PCHMG));
            PetscCallThrow(PCSetReusePreconditioner(pc, PETSC_FALSE));

            PetscCallThrow(PCHMGSetInnerPCType(pc, PCGAMG));
            PetscCallThrow(PCHMGSetReuseInterpolation(pc, PETSC_FALSE));
            PetscCallThrow(PCHMGSetUseSubspaceCoarsening(pc, PETSC_TRUE));
            PetscCallThrow(PCHMGUseMatMAIJ(pc, PETSC_FALSE));
            PetscCallThrow(PCHMGSetCoarseningComponent(pc, 0));
        }
    }


    //PetscCallThrow(KSPSetOperators(ksp, K, K));
    //PetscCallThrow(KSPSetUp(ksp));
    PetscCallThrow(KSPSolve(ksp, f, dx));


    firstSolve = false;
}