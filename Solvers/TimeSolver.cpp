
#include "TimeSolver.h"
#include "../InputsOutputs/Restarter.h"

/// @brief Initialization based on input file
/// @param inputs Input file object
/// @param physics_in physics object to link time sovler to
/// @param NLSolver_in nonlinear solver to use within time solver
/// @param Restart_in Object to save restartable data
TimeSolver::TimeSolver(inputData& inputs, Physics& physics_in, NonLinSolver& NLSolver_in, Restarter& Restart_in){
    Setup(inputs, physics_in, NLSolver_in, Restart_in);

    step = 0;
    physics->time = t;
    physics->timeScheme->Set_dt(dt);
    physics->MaxTime = tmax;
}

/// @brief Restart from a previous simulation
/// @param inputs Input file object
/// @param physics_in physics object to link time sovler to
/// @param NLSolver_in nonlinear solver to use within time solver
/// @param Restart_in Object to save restartable data
/// @param data Restart data
TimeSolver::TimeSolver(inputData& inputs, Physics& physics_in, NonLinSolver& NLSolver_in, Restarter& Restart_in, SaveDataFile& data){
    Setup(inputs, physics_in, NLSolver_in, Restart_in);

    data.Load("Step", step);
    data.Load("t", t);
    data.Load("dt",dt);
    physics->time = t;
    physics->timeScheme->Set_dt(dt);
    physics->MaxTime = tmax;
}

/// @brief Setup required independent of start method
/// @param inputs input file
/// @param physics_in physics object to link time sovler to
/// @param NLSolver_in nonlinear solver to use within time solver
/// @param Restart_in Object to save restartable data
void TimeSolver::Setup(inputData& inputs, Physics& physics_in, NonLinSolver& NLSolver_in, Restarter& Restart_in){
    PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));

    physics = &physics_in;
    NRSolver = &NLSolver_in;
    RestartSaver = &Restart_in;

    // initial time
    t=0.0;
    bool hast0 = inputs.GetOptional(t, {"TimeSolver","t0"});

    //time icnrements

    UseDtSteps = false;
    UseDtSteps = inputs.GetOptional(dtSteps, {"TimeSolver","dtSteps"});
    if (UseDtSteps){
        inputs.GetRequired(tSteps, {"TimeSolver","tSteps"});
    } else {
        inputs.GetRequired(dt, {"TimeSolver","dt"});
        if (hast0){
            inputs.GetRequired(dt0,{"TimeSolver","dt0"});
        } else {
            dt0 = 1.0;
        }

        bool hastGrow = inputs.GetOptional(dtGrow, {"TimeSolver","dtGrow"});
        if (hastGrow){
            inputs.GetRequired(dtMax,{"TimeSolver","dtMax"});
        } else {
            dtGrow=1.0;
            dtMax = dt;
        } 
    }

    DynamicDT = false;
    inputs.GetOptional(DynamicDT, {"TimeSolver","DynamicDT"});

    inputs.GetRequired(tmax, {"TimeSolver","tmax"});
    inputs.GetRequired(outputN, {"TimeSolver","outputN"});

    RefineByData = inputs.GetOptional(TimeRefineDataName, {"TimeSolver","TimeDataRefine"});
    if (RefineByData){
        inputs.GetRequired(TimeRefineDataLim, {"TimeSolver","TimeDataRefineFactor"});

        // find timedata Type
        uint TimeRefineData;
        TimeRefineData = physics->TimeDataTypes.size();
        for (size_t i = 0; i < physics->TimeDataTypes.size(); i++){
            if (TimeRefineDataName == physics->TimeDataTypes[i]){
                TimeRefineData = i;
            }
        }
        if (TimeRefineData == physics->TimeDataTypes.size()){
            throw std::invalid_argument("TimeRefineData "+TimeRefineDataName+" not found in physics->TimeDataTypes");
        }
    }

    Saver = new OutputDataSaver(physics_in, inputs);
    Saver->saveResults(0, t);
}

/// @brief Save to restart file
/// @param data Data object to save into
void TimeSolver::Save(SaveDataFile& data){
    data.Save("Step", step);
    data.Save("t", t);
    data.Save("dt",dt);
}

TimeSolver::~TimeSolver(){
    delete Saver;
}

/// @brief Solves the time-dependent problem
void TimeSolver::Solve(){
    RestartSaver->Save(0);
    dtScale = 1.0;

    while (t<physics->MaxTime){
        Logs.PrintSingle("t="+std::to_string(t)+" (timestep "+std::to_string(step)+"):\n",0);
        physics->time = t;

        Setdt();

        bool StepSolved, hasConverged;

        StepSolved = false;
        while (StepSolved==false){
            try{
                hasConverged = NRSolver->SolveStep();
            } catch(const std::exception& e){
                hasConverged = false;
                if (DynamicDT && dtScale>1.0e-9){
                    std::cout << e.what() << '\n';
                } else {
                    throw;
                }
            }
            bool RefineStep = false;
            if (hasConverged==false && DynamicDT && dtScale>1.0e-9){
                Logs.PrintSingle("\n\nStep did not converge, reducing time step\n",1);
                RefineStep = true;
            } else if (RefineByData && DynamicDT && dtScale>1.0e-9){
                double tData = physics->GetTimeData(TimeRefineDataName);
                if (tData>TimeRefineDataLim && t>-1.0e-3){
                    std::ostringstream Msg; Msg << "TimeRefineData exceeded limit ("<< TimeRefineDataName << ": " << tData << ">" << TimeRefineDataLim << "), reducing time step\n";
                    Logs.PrintSingle(Msg,1);
                    RefineStep = true;
                    hasConverged = false;
                } else {
                    StepSolved = true;
                }
            } else {
                StepSolved = true;
            }
            if (RefineStep){
                MPI_Barrier(PETSC_COMM_WORLD);
                dtScale = dtScale/2.0;
                Setdt();
                physics->ResetState();
                std::ostringstream Msg2; Msg2 << "Reduced dt to " << physics->timeScheme->dt << "\n\n";
                Logs.PrintSingle(Msg2,1);
            }
        }
        

        physics->Commit(CommitTypes::TIMEDEP_COMMIT_TYPE);
        t += physics->timeScheme->dt;
        step += 1;

        if (dtScale<1.0){
            dtScale *= 1.1;
            dtScale = std::min(1.0, dtScale);
        }

        if ( step%outputN == 0){
            Saver->saveResults(step, t);
        }
        RestartSaver->Save(step);
    }
}

void TimeSolver::Setdt(){
    if (UseDtSteps){
        for (size_t i = 0; i < dtSteps.size(); i++){
            if (t>=tSteps[i]){
                physics->timeScheme->Set_dt(dtScale*dtSteps[i]);
            }
        }
    } else {
        if (t<0.0){
            physics->timeScheme->Set_dt(dtScale*dt0);
        } else {
            if (t>0 && dtGrow>1.0){
                dt = std::min(dtMax, dt*dtGrow);
            }
            physics->timeScheme->Set_dt(dtScale*dt);
        }
    }
}