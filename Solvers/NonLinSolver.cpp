#include <chrono>

#include "NonLinSolver.h"

/// @brief Creates a nonlinear solver based on an input file
/// @param inputs Input file
/// @param physics_in Reference to hpysics object to solve for
NonLinSolver::NonLinSolver(inputData& inputs, Physics& physics_in){
    Setup(inputs, physics_in);
}

/// @brief Restarts a nonlinear solver based on an input and restart file
/// @param inputs Input file
/// @param physics_in Reference to hpysics object to solve for
/// @param data Restart data file
NonLinSolver::NonLinSolver(inputData& inputs, Physics& physics_in, SaveDataFile& data){
    Setup(inputs, physics_in);
}

/// @brief Save solver to restart file
/// @param data Restart data file to save into
void NonLinSolver::Save(SaveDataFile& data){

}

/// @brief Performs setup required to be done before use
/// @param inputs iput data file
/// @param physics_in Reference to physics object
void NonLinSolver::Setup(inputData& inputs, Physics& physics_in){
    PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));

    physics = &physics_in;

    //Load data
    inputs.GetRequired(maxIt, {"nonlinSolver","max_it"});
    inputs.GetRequired(tiny, {"nonlinSolver","tiny"});
    inputs.GetRequired(max_outer_it, {"nonlinSolver","max_outer_it"});


    std::string LinSolveType;  inputs.GetRequired(LinSolveType,{"nonlinSolver","LinSolver","Type"});
    nSteps = physics->dofspace->maxSteps;
    inputs.GetRequired(convCrit, {"nonlinSolver","convCrit"});

    bool hasLineSearchOption = inputs.GetOptional(lineSearch, {"nonlinSolver","linesearch"});
    if (hasLineSearchOption==false){
        lineSearch = false;
    }
    if (lineSearch){
        inputs.GetRequired(lineSearchLimits,{"nonlinSolver","linesearchLims"});
        if (lineSearchLimits.size()!=2){
            throw std::invalid_argument("nonlinsolver.linesearchLims should specify both an upper and lower limit");
        }
    }

    //resize vectors
    dx.resize(nSteps);
    Err0.resize(nSteps);
    ErrIt0.resize(nSteps);
    Err.resize(nSteps);

    for (size_t i = 0; i < nSteps; i++){
        VecCreateMPI(MPI_COMM_WORLD, physics->cons->ConstrainedSizeLocal[i], PETSC_DECIDE, &dx[i]);
        VecZeroEntries(dx[i]);

        Err0[i].resize(3);
        ErrIt0[i].resize(3);
        Err[i].resize(3);
    }

    //create linear solver
    LinSolvers.resize(nSteps);
    for (size_t i = 0; i < nSteps; i++){
        LinSolvers[i] = CreateLinSolver(LinSolveType);
        LinSolvers[i]->init(inputs, physics->K[i]);
    }

    //BFGS
    BFGS.resize(nSteps);
    for (size_t i = 0; i < nSteps; i++){
        BFGS[i] = false;
    }
    inputs.GetOptional(BFGS, {"nonlinSolver","BFGS"});
    if (std::find(begin(BFGS), end(BFGS), true) != end(BFGS)){ // At least one step uses BFGS
        dX.resize(nSteps);
        dr.resize(nSteps);
        bfgs_a.resize(nSteps);
        bfgs_b.resize(nSteps);
        bfgs_rho.resize(nSteps);
        q.resize(nSteps);
    }
    for (size_t i = 0; i < nSteps; i++){
        if (BFGS[i]){
            BFGSSetup(inputs, i);
        }
    }

    Saver = new OutputDataSaver(physics_in, inputs);
}

NonLinSolver::~NonLinSolver(){
    for (size_t i = 0; i < nSteps; i++){
        delete LinSolvers[i];
    }
}

/// @brief Solves a single time step
bool NonLinSolver::SolveStep(){
    bool convergedGlobal = false, stopGlobal = false;
    std::vector<bool> ConvergedSteps(nSteps, false);

    for (size_t i = 0; i < nSteps; i++){
        Err0[i][0] = -1.0; Err0[i][1] = -1.0; Err0[i][2] = -1.0;
    }
    
    size_t outerIt = 0;
    while (convergedGlobal == false && stopGlobal == false){
        for (size_t i = 0; i < nSteps; i++){
            ErrIt0[i][0] = -1.0; ErrIt0[i][1] = -1.0; ErrIt0[i][2] = -1.0;
        }

        Logs.PrintSingle("Loop "+std::to_string(outerIt)+"\n",1);
        for (size_t step = 0; step < nSteps; step++){
            Logs.PrintSingle("Step " + std::to_string(step) + "\n", 1);
            if (BFGS[step]){
                ConvergedSteps[step] = SolveStep_BFGS(step);
            } else {
                ConvergedSteps[step] = SolveStep_Standard(step);
            }
        }
        outerIt += 1;

        //check global convergence
        convergedGlobal = true;
        for (size_t i = 0; i < nSteps; i++){
            if (ConvergedSteps[i] == false){
                convergedGlobal = false;
            } else if (max_outer_it>1){
                for (size_t j = 0; j < convCrit.size(); j++){
                    if (ErrIt0[i][j] > convCrit[j] && ErrIt0[i][j]*Err0[i][j]>tiny[j]){
                        convergedGlobal = false;
                    }
                } 
            }
        }

        if (outerIt>=max_outer_it){
            stopGlobal = true;
        }
    }
    return convergedGlobal;
}

bool NonLinSolver::SolveStep_Standard(size_t step){
    std::chrono::high_resolution_clock::time_point start, end;

    size_t it = 0;
    bool converged = false, stop = false;
    physics->Assemble(step);
    physics->Constraint(step);
    while (converged == false && stop == false){
        Logs.PrintSingle("Iteration " + std::to_string(it) + "\n", 1);

        Logs.PrintSingle("Solving\n", 2);
        start = std::chrono::high_resolution_clock::now();
        LinSolvers[step]->solve(physics->KCon[step], physics->fCon[step], dx[step]);
        end = std::chrono::high_resolution_clock::now();
        Logs.PrintSingle("Solver Time: " + std::to_string(std::chrono::duration<double, std::milli>(end - start).count() / 1000.0) + "s\n", 3);

        if (lineSearch && (it > 0)){
            HOLineSearch(step);
        }
        else{
            physics->UpdateState(step, dx[step]);
        }
        physics->Assemble(step);
        physics->Constraint(step);

        converged = CheckConvergence(step);

        // Saver->saveResults(0,0.0);

        #ifdef LOCALENVIRONMENT
            physics->Vis_Plots->Plot("it", step);
        #endif

        it += 1;
        if (it > maxIt){
            stop = true;
        }
    }
    return converged;
}

/// @brief Check convergence criteria for the current staggered solver step
/// @param step Staggered solver step
/// @return convergence criteria fulfilled
bool NonLinSolver::CheckConvergence(size_t step){
    bool converged = true;

    Vec e;
    PetscInt nDofs; 
    VecDuplicate(dx[step], &e);
    VecGetSize(e, &nDofs);

    VecNorm(dx[step], NORM_1, &Err[step][0]); //displacement increment norm
    VecNorm(physics->fCon[step], NORM_1, &Err[step][1]); //force mismatch norm
    VecPointwiseMult(e, physics->fCon[step], dx[step]);
    VecNorm(e, NORM_1, &Err[step][2]); //energy mismatch
    VecDestroy(&e);

    for (size_t i = 0; i < 3; i++){
        Err[step][i] = Err[step][i]/nDofs;
    }
    

    if (Err0[step][0] < 0){
        Err0[step][0] = Err[step][0];
        Err0[step][1] = Err[step][1];
        Err0[step][2] = Err[step][2];
        if (Err0[step][0]<tiny[0]) Err0[step][0] = tiny[0];
        if (Err0[step][1]<tiny[1]) Err0[step][1] = tiny[1];
        if (Err0[step][2]<tiny[2]) Err0[step][2] = tiny[2];
    }
    for (size_t i = 0; i < 3; i++){
        Err[step][i] /= Err0[step][i];

        if (Err[step][i] > convCrit[i] && Err[step][i]*Err0[step][i]>tiny[i]){
            converged = false;
        }
    }
    if (ErrIt0[step][0] < 0){
        ErrIt0[step][0] = Err[step][0];
        ErrIt0[step][1] = Err[step][1];
        ErrIt0[step][2] = Err[step][2];
    }

    std::ostringstream Msg; Msg << "Errors: ";
    for (size_t i = 0; i < 3; i++){
        Msg << Err[step][i] << "  ("  << Err[step][i]*Err0[step][i] << "), ";
    }
    Msg << "\n";

    bool hasNaNs = false;
    if (isnan(Err[step][0]) || Err[step][0] > 1e6){
        Msg << "Displacement increment norm NaN, exiting\n";
        hasNaNs = true;
    }
    if (isnan(Err[step][1] || Err[step][1] > 1e6)){
        Msg << "Force norm NaN, exiting\n";
        hasNaNs = true;
    }
    if (isnan(Err[step][2] || Err[step][2] > 1e6)){
        Msg << "Energy norm NaN, exiting\n";
        hasNaNs = true;
    }

    Logs.PrintSingle(Msg,2);
    if (hasNaNs){
        //Saver->saveResults(0,0.0);
        PetscCallThrow(PetscBarrier(NULL));
        throw std::invalid_argument("NaN in Force/displacement/Energy residuals \n");
    }
    
    return converged;
}

/// @brief Perform a linear line-search
/// @param step Staggered solver step
void NonLinSolver::LinLineSearch(size_t step){
    std::ostringstream Msg; Msg << "Linesearch: ";
    double e0, e1;
    double L_Search;
    Vec dxCopy;
    VecDuplicate(dx[step], &dxCopy);
    VecCopy(dx[step], dxCopy);

    //error before step
    VecDot(physics->fCon[step], dxCopy, &e0);
    Msg << e0 << " -> ";

    //error after step
    physics->UpdateState(step, dxCopy);
    physics->AssembleFOnly(step);
    physics->ConstraintFOnly(step);
    VecDot(physics->fCon[step], dxCopy, &e1);
    Msg << e1 << ", Factor:";

    //perform line-search
    if ((e1-e0)==0.0){
        L_Search = 1.0;
    } else {
        L_Search = -e0/(e1-e0);
        if (L_Search < lineSearchLimits[0]) L_Search = lineSearchLimits[0];
        if (L_Search > lineSearchLimits[1]) L_Search = lineSearchLimits[1];
    }
    VecScale(dxCopy, -(1.0-L_Search));
    VecScale(dx[step], L_Search);
    physics->UpdateState(step, dxCopy);
    Msg << L_Search << "\n";
    Logs.PrintSingle(Msg,2);

    VecDestroy(&dxCopy);
}

void NonLinSolver::HOLineSearch(size_t step){
    uint maxLS = 10;
    bool found = false;
    double e0;
    double L_Search_corr, eMin;
    std::vector<double> e(maxLS+2), L_Search(maxLS+2);
    Vec dxCopy; VecDuplicate(dx[step], &dxCopy);
    VecDot(physics->fCon[step], dx[step], &e0);

    size_t it = 0;
    L_Search[0] = 1.0;
    std::ostringstream Msg; Msg << "Linesearch: e0=" << e0 << "\n";
    Logs.PrintSingle(Msg,2);
    while (found == false){
        VecCopy(dx[step], dxCopy);
        VecScale(dxCopy, L_Search[it]);
        physics->UpdateState(step, dxCopy);

        physics->AssembleFOnly(step);
        physics->ConstraintFOnly(step);
        VecDot(physics->fCon[step], dxCopy, &e[it]);

        std::ostringstream Msg; Msg << "Linesearch: ";
        Msg << "it " << it << ": L=" << L_Search[it] << "(e=" << e[it] << ")\n";
        Logs.PrintSingle(Msg,2);

        if ((std::abs(e[it]) < std::abs(e0) && sgn(e[it]) == sgn(e0)) || e[it] == 0.0){
            found = true;
            L_Search_corr = L_Search[it];
            eMin = e[it];
        } else {
            if (sgn(e[it]) != sgn(e0)){
                L_Search[it+1] = -e0/(e[it]-e0)*L_Search[it];
                L_Search[it+1] = std::max(lineSearchLimits[0], std::min(lineSearchLimits[1], L_Search[it+1]));
            } else {
                L_Search[it+1] = 0.5*L_Search[it];
            }

            if (std::abs(L_Search[it+1]-L_Search[it]) <= lineSearchLimits[0]){
                found = true;
                L_Search_corr = L_Search[it];
                eMin = e[it];
            }
        }

        VecCopy(dx[step], dxCopy);
        VecScale(dxCopy, -L_Search[it]);
        physics->UpdateState(step, dxCopy);

        it += 1;
        if (it >= maxLS){
            found = true;
            eMin = e0;
            L_Search_corr = lineSearchLimits[0];
            for (size_t i = 0; i < it-1; i++){
                if (e[i] < eMin){
                    eMin = e[i];
                    L_Search_corr = L_Search[i];
                }
            }
        }
    }

    VecScale(dx[step], L_Search_corr);
    physics->UpdateState(step, dx[step]);

    std::ostringstream Msg2; Msg2 << "Linesearch final: ";
    Msg2 << L_Search_corr << "(" << e0 << "->" << eMin << ")\n";
    Logs.PrintSingle(Msg2,2);

    VecDestroy(&dxCopy);
}

/// @brief Set up BFGS solver for a given staggered step
/// @param inputs Input data file
/// @param step Staggered step number
void NonLinSolver::BFGSSetup(inputData& inputs, size_t step){
    BFGSReset = maxIt;
    inputs.GetOptional(BFGSReset, {"nonlinSolver","BFGS_Reset"});

    bfgs_a[step].resize(BFGSReset);
    bfgs_b[step].resize(BFGSReset);
    bfgs_rho[step].resize(BFGSReset);
    dX[step].resize(BFGSReset);
    dr[step].resize(BFGSReset);

    VecCreateMPI(MPI_COMM_WORLD, physics->cons->ConstrainedSizeLocal[step], PETSC_DECIDE, &q[step]);
    for (size_t i = 0; i < BFGSReset; i++){
        VecCreateMPI(MPI_COMM_WORLD, physics->cons->ConstrainedSizeLocal[step], PETSC_DECIDE, &dX[step][i]);
        VecCreateMPI(MPI_COMM_WORLD, physics->cons->ConstrainedSizeLocal[step], PETSC_DECIDE, &dr[step][i]);
    }
}

/// @brief Perform a BFGS step for a given staggered step
/// @param step Staggered step number
bool NonLinSolver::SolveStep_BFGS(size_t step){
    std::chrono::high_resolution_clock::time_point start, end;

    size_t it = 0;
    size_t bfgs_it = 0;
    size_t bfgs_failed = 0;
    bool converged = false, stop = false;
    physics->Assemble(step);
    physics->Constraint(step);
    while (converged == false && stop == false){
        Logs.PrintSingle("Iteration " + std::to_string(it) + "(BFGS it "+ std::to_string(bfgs_it) + ")\n", 1);

        Logs.PrintSingle("Solving\n", 2);
        start = std::chrono::high_resolution_clock::now();

        if (bfgs_it==0){
            LinSolvers[step]->solve(physics->KCon[step], physics->fCon[step], dx[step]);
        } else {
            VecCopy(physics->fCon[step], q[step]);
            VecScale(q[step], -1.0);

            double ai, bi;
            for (int i = bfgs_it-1; i >= 0; i--){
                VecDot(q[step], dX[step][i], &ai);
                bfgs_a[step][i] = bfgs_rho[step][i]*ai;
                
                VecAXPY(q[step], -bfgs_a[step][i], dr[step][i]);
            }
            LinSolvers[step]->solve(physics->KCon[step], q[step], dx[step]);
            for (size_t i = 0; i <= bfgs_it-1; i++){
                VecDot(dx[step], dr[step][i], &bi);
                bfgs_b[step][i] = bfgs_rho[step][i]*bi;
                VecAXPY(dx[step], bfgs_a[step][i] - bfgs_b[step][i], dX[step][i]);
            }
            VecScale(dx[step], -1.0);
        }
        
        end = std::chrono::high_resolution_clock::now();
        Logs.PrintSingle("Solver Time: " + std::to_string(std::chrono::duration<double, std::milli>(end - start).count() / 1000.0) + "s\n", 3);


        VecCopy(physics->fCon[step], dr[step][bfgs_it]);

        if (lineSearch && (it > 0)){
            HOLineSearch(step);
        }else{
            physics->UpdateState(step, dx[step]);
        }

        if (bfgs_it < BFGSReset-1 && bfgs_failed < 5){
            physics->AssembleFOnly(step);
            physics->ConstraintFOnly(step);

            VecCopy(dx[step], dX[step][bfgs_it]);
            VecAXPY(dr[step][bfgs_it], -1.0, physics->fCon[step]);

            double rhoi;
            VecDot(dX[step][bfgs_it], dr[step][bfgs_it], &rhoi);
            if (rhoi<=0.0){
                //increment not added
                Logs.PrintSingle("BFGS increment not added\n", 2);
                bfgs_failed += 1;
            } else {
                //increment added
                bfgs_rho[step][bfgs_it] = 1.0/rhoi;
                bfgs_it += 1;
            }
        } else {
            bfgs_it = 0;
            bfgs_failed = 0;
            physics->Assemble(step);
            physics->Constraint(step);
        }

        converged = CheckConvergence(step);

        #ifdef LOCALENVIRONMENT
            physics->Vis_Plots->Plot("it", step);
        #endif

        it += 1;
        if (it > maxIt){
            stop = true;
        }
    }

    return converged;
}