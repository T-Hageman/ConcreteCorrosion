#include "physics.h"
#include "../utility/utility.h"
#include "../Initializations/Initializations.h"

#include <petscviewerhdf5.h>

#include <petsc.h>
#include <iostream>
#include <string>
#include <chrono>

/// @brief Initialize physics based on an input file
/// @param inputs Input file
/// @param inMesh reference to mesh
Physics::Physics(inputData& inputs, Mesh& inMesh) {
    PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));

    //initialize objects
    mesh = &inMesh;
    dofspace = new DofSpace(inputs);
    cons = new Constrainer(dofspace);
    timeScheme = new TimeScheme(inputs); timeScheme->Set_dt(1.0);
    time = 0;

    //initialize physics models
    inputs.GetRequired(ModelNames, {"Models","Names"});
    nModels = ModelNames.size();
    
    ModelTypeNames.resize(nModels);
    Models.resize(nModels);
    for (size_t i = 0; i < nModels; i++){
        Logs.PrintSingle("Initializing:"+ModelNames[i]+"\n",2);
        inputs.GetRequired(ModelTypeNames[i], {"Models", ModelNames[i], "Name"});
        Models[i] = CreateModel(*this, ModelTypeNames[i], ModelNames[i]);

        Models[i]->init(inputs);
        Logs.PrintSingle("Finished initializing:"+ModelNames[i]+"\n",2);
    }
    dofspace->SyncDofs(mesh);
    dofspace->ExportDofSpace();

    //set up state vectors
    InitStateVecs();
    InitKf();
    Initialize(inputs, *this);
    for (size_t i = 0; i < dofspace->maxSteps; i++){
        StateVectorsOld[i].Set(StateVectors[i]);
        dStateVectorsOld[i].Set(dStateVectors[i]);
        ddStateVectorsOld[i].Set(ddStateVectors[i]);

        StateVectorsOld[i].AssemblyStart();
        StateVectorsOld[i].AssemblyEnd();
        StateVectorsOld[i].SyncStart(INSERT_VALUES);
        StateVectorsOld[i].SyncEnd(INSERT_VALUES);

        dStateVectorsOld[i].AssemblyStart();
        dStateVectorsOld[i].AssemblyEnd();
        dStateVectorsOld[i].SyncStart(INSERT_VALUES);
        dStateVectorsOld[i].SyncEnd(INSERT_VALUES);

        ddStateVectorsOld[i].AssemblyStart();
        ddStateVectorsOld[i].AssemblyEnd();
        ddStateVectorsOld[i].SyncStart(INSERT_VALUES);
        ddStateVectorsOld[i].SyncEnd(INSERT_VALUES);
    }

    //test assemble to see if all correct
    Logs.PrintSingle("Setting up matrices for first time, this takes a moment\n",2);
    for (size_t i = 0; i < dofspace->maxSteps; i++){
        Logs.PrintSingle("Assembling step "+std::to_string(i)+"\n",2);
        Assemble(i);
        MPI_Barrier(PETSC_COMM_WORLD);
        Logs.PrintSingle("Adding constraints step "+std::to_string(i)+"\n",2);
        Constraint(i);
        MPI_Barrier(PETSC_COMM_WORLD);
        Logs.PrintSingle("Finished assembling step "+std::to_string(i)+"\n",2);
    }
    BoundaryConsToOld();

    Init_Timeseries();
    Logs.PrintSingle("Finished setting up matrices\n",2);

    #ifdef LOCALENVIRONMENT
        Vis_Plots = new Visualisation(inputs, this);
        Vis_Plots->Plot("step", 0);
    #endif
}

/// @brief Restart physics from previosuly run simulation
/// @param inputs Input data file
/// @param inMesh reference to mesh
/// @param data Restart data file
Physics::Physics(inputData& inputs, Mesh& inMesh, SaveDataFile& data){
    PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
    std::string OldPrefix = data.GetPrefix();

    double dt;

    mesh = &inMesh;
    data.Load("time", time);
    data.Load("dt", dt);
    data.SetPrefix(OldPrefix+"/DofSpace");
    dofspace = new DofSpace(inputs);
    dofspace->load(data);
    cons = new Constrainer(dofspace);
    timeScheme = new TimeScheme(inputs);

    timeScheme->Set_dt(dt);

    inputs.GetRequired(ModelNames, {"Models","Names"});
    nModels = ModelNames.size();

    ModelNames.resize(nModels); ModelTypeNames.resize(nModels);
    Models.resize(nModels);

    for (size_t i = 0; i < nModels; i++){
        inputs.GetRequired(ModelTypeNames[i], {"Models",ModelNames[i],"Name"});
        Models[i] = CreateModel(*this, ModelTypeNames[i], ModelNames[i]);

        data.SetPrefix(OldPrefix+"/"+ModelNames[i]);
        Models[i]->load(inputs, data);
    }
    data.SetPrefix(OldPrefix);

    size_t nSteps = dofspace->maxSteps;
    StateVectors.resize(nSteps);
    StateVectorsOld.resize(nSteps);
    dStateVectors.resize(nSteps);
    ddStateVectors.resize(nSteps);
    dStateVectorsOld.resize(nSteps);
    ddStateVectorsOld.resize(nSteps);
    for (size_t i = 0; i < nSteps; i++){
        data.Load("StateVector_"+std::to_string(i), StateVectors[i]);
        data.Load("StateVectorsOld_"+std::to_string(i), StateVectorsOld[i]);
        data.Load("dStateVectors_"+std::to_string(i), dStateVectors[i]);
        data.Load("dStateVectorsOld_"+std::to_string(i), dStateVectorsOld[i]);
        data.Load("ddStateVectors_"+std::to_string(i), ddStateVectors[i]);
        data.Load("ddStateVectorsOld_"+std::to_string(i), ddStateVectorsOld[i]);
    }
    InitKf();
    for (size_t i = 0; i < dofspace->maxSteps; i++){
        Assemble(i);
        Constraint(i);
    }

    Init_Timeseries();
    data.Load("TimeData", TimeData);

    #ifdef LOCALENVIRONMENT
        Vis_Plots = new Visualisation(inputs, this);
        Vis_Plots->Plot("step", 0);
    #endif
}

/// @brief Save data required for restarts
/// @param data Restart data object into which the physics model is saved
void Physics::Save(SaveDataFile& data){
    std::string OldPrefix = data.GetPrefix();
    data.Save("time", time);
    data.Save("dt", timeScheme->dt);

    for (size_t i = 0; i < nModels; i++){
        data.SetPrefix(OldPrefix+"/"+ModelNames[i]);
        Models[i]->save(data);
    }

    data.SetPrefix(OldPrefix+"/DofSpace");
    dofspace->save(data);

    data.SetPrefix(OldPrefix);
    size_t nSteps = dofspace->maxSteps;
    for (size_t i = 0; i < nSteps; i++){
        data.Save("StateVector_"+std::to_string(i), StateVectors[i]);
        data.Save("StateVectorsOld_"+std::to_string(i), StateVectorsOld[i]);
        data.Save("dStateVectors_"+std::to_string(i), dStateVectors[i]);
        data.Save("dStateVectorsOld_"+std::to_string(i), dStateVectorsOld[i]);
        data.Save("ddStateVectors_"+std::to_string(i), ddStateVectors[i]);
        data.Save("ddStateVectorsOld_"+std::to_string(i), ddStateVectorsOld[i]);
    }

    data.Save("TimeData", TimeData);
}

Physics::~Physics() {
    for (size_t i = 0; i < nModels; i++){
        delete Models[i];
    }
    
    delete dofspace;
    delete cons;
    delete timeScheme;
}

/// @brief Initialize, resize and zero state vectors
void Physics::InitStateVecs(){
    size_t nSteps = dofspace->maxSteps;
    StateVectors.resize(nSteps);
    StateVectorsOld.resize(nSteps);
    dStateVectors.resize(nSteps);
    ddStateVectors.resize(nSteps);
    dStateVectorsOld.resize(nSteps);
    ddStateVectorsOld.resize(nSteps);
    for (size_t dStep = 0; dStep < nSteps; dStep++){
        size_t localDofCount = dofspace->LocalDofNumRange[dStep][1]-dofspace->LocalDofNumRange[dStep][0];
        size_t globalDofCount = dofspace->TotalDofRange[dStep];
        std::vector<PetscInt> GhostDofs = dofspace->getGhostDofNumbers(dStep);

        StateVectors[dStep].SetSize(localDofCount, globalDofCount, GhostDofs);
        StateVectors[dStep].Zero();
        StateVectorsOld[dStep].SetSize(localDofCount, globalDofCount, GhostDofs);
        StateVectorsOld[dStep].Zero();

        dStateVectors[dStep].SetSize(localDofCount, globalDofCount, GhostDofs);
        dStateVectors[dStep].Zero();
        dStateVectorsOld[dStep].SetSize(localDofCount, globalDofCount, GhostDofs);
        dStateVectorsOld[dStep].Zero();

        ddStateVectors[dStep].SetSize(localDofCount, globalDofCount, GhostDofs);
        ddStateVectors[dStep].Zero();
        ddStateVectorsOld[dStep].SetSize(localDofCount, globalDofCount, GhostDofs);
        ddStateVectorsOld[dStep].Zero();
    }
}

/// @brief Initialize force vector and stiffness matrx for all staggered steps
void Physics::InitKf(){
    size_t nSteps = dofspace->maxSteps;
    f.resize(nSteps); fCon.resize(nSteps);
    K.resize(nSteps); KCon.resize(nSteps);
    CForceVec.resize(nSteps); CForceVec2.resize(nSteps);
    dofVersionUsed.resize(nSteps);
    for (size_t dStep = 0; dStep < nSteps; dStep++){
        size_t localDofCount = dofspace->LocalDofNumRange[dStep][1]-dofspace->LocalDofNumRange[dStep][0];
        size_t globalDofCount = dofspace->TotalDofRange[dStep];
        
        PetscCallThrow(VecCreateMPI(PETSC_COMM_WORLD, localDofCount, globalDofCount, &f[dStep]));
        PetscCallThrow(VecCreateMPI(PETSC_COMM_WORLD, localDofCount, globalDofCount, &CForceVec[dStep]));
        PetscCallThrow(VecCreateMPI(PETSC_COMM_WORLD, localDofCount, globalDofCount, &CForceVec2[dStep]));
        PetscCallThrow(MatCreateAIJ(PETSC_COMM_WORLD, localDofCount, localDofCount, globalDofCount, globalDofCount, nAlloc, NULL, nAlloc, NULL, &K[dStep]));

        MatSetOption(K[dStep], MAT_ROW_ORIENTED, PETSC_FALSE);
        MatSetOption(K[dStep], MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
        
        PetscCallThrow(VecZeroEntries(f[dStep]));
        PetscCallThrow(MatZeroEntries(K[dStep]));

        dofVersionUsed[dStep] = 0;
    }
}

/// @brief Assemble the force vector and tangential matrix
/// @param step current staggered step number
void Physics::Assemble(size_t step){
    Logs.PrintSingle("Starting matrix assembly:\n",2);

    // sync state vectors across cores and set force and stiffness to zero
    StateVectors[step].SyncStart(INSERT_VALUES);
    PetscCallThrow(VecZeroEntries(f[step]));
    PetscCallThrow(MatZeroEntries(K[step]));
    cons->SetZero(step);
    StateVectors[step].SyncEnd(INSERT_VALUES);

    //loop over models to calculate tangential matrix
    std::chrono::high_resolution_clock::time_point start, end;
    std::vector<double> Timings(nModels), MaxTimings(nModels), MinTimings(nModels);
    std::vector<MPI_Request> MPI_Reqs(2);
    for (size_t i = 0; i < nModels; i++){
        Logs.PrintSingle(Models[i]->MyName+":\n",3);
        start = std::chrono::high_resolution_clock::now();
        Models[i]->Assemble(K[step], f[step], cons, step);
        end = std::chrono::high_resolution_clock::now();
        Timings[i] = std::chrono::duration<double, std::milli>( end - start ).count();
        //MPI_Ireduce(&Timings[i], &MaxTimings[i], 1,MPI_DOUBLE, MPI_MAX, 0, PETSC_COMM_WORLD, &MPI_Reqs[0]);
        //MPI_Ireduce(&Timings[i], &MinTimings[i], 1,MPI_DOUBLE, MPI_MIN, 0, PETSC_COMM_WORLD, &MPI_Reqs[1]);
        if (rank==0){
            //std::vector<MPI_Status> myStatus(MPI_Reqs.size());
            //MPI_Waitall(MPI_Reqs.size(), MPI_Reqs.data(), myStatus.data());
            //PrintSingle("\t\t\t\tRuntime: min "+std::to_string(MinTimings[i]/1000.0)+"s  max "+std::to_string(MaxTimings[i]/1000.0)+"s\n");
            Logs.PrintSingle("\tRuntime: min "+std::to_string(Timings[i]/1000.0)+"s  max "+std::to_string(Timings[i]/1000.0)+"s\n",3);
        }
    }

    //final asembly and cross-core communications
    cons->Assemble(step, StateVectors[step]);
    PetscCallThrow(VecAssemblyBegin(f[step]));
    PetscCallThrow(MatAssemblyBegin(K[step], MAT_FINAL_ASSEMBLY));
    PetscCallThrow(VecAssemblyEnd(f[step]));
    PetscCallThrow(MatAssemblyEnd(K[step], MAT_FINAL_ASSEMBLY));
    cons->AssembleEnd(step);
    PetscCallThrow(PetscBarrier(NULL));

    //debugging function, export matrix to file for viewing within matlab
    if (exportMatrices){
        PetscViewer MatViewer;
        std::string FName = "mat"+std::to_string(step)+".bin";
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, FName.c_str(), FILE_MODE_WRITE, &MatViewer);
        PetscCallThrow(MatView(K[step], MatViewer));
        PetscCallThrow(PetscViewerDestroy(&MatViewer));
    }
}

/// @brief Assemble the force vector only
/// @param step current staggered step number
void Physics::AssembleFOnly(size_t step){
    Assemble(step);
    //Still to update to only do F
}

/// @brief Apply constraints to the force vector for the current staggered step
/// @param step staggered step number
void Physics::ConstraintFOnly(size_t step){
    if (dofVersionUsed[step] != cons->dofVersion[step]){
        PetscCallThrow(VecCreateMPI(PETSC_COMM_WORLD, cons->ConstrainedSizeLocal[step], PETSC_DECIDE, &fCon[step]));
    }
    PetscCallThrow(MatMult(cons->ConMats[step], cons->ConValVec[step], CForceVec[step]));   //  F_c = C*C_v
    PetscCallThrow(VecScale(CForceVec[step], -1.0));
    PetscCallThrow(MatMultAdd(K[step], CForceVec[step], f[step], CForceVec2[step]));        //  F_c'= f-K*C*C_v
    PetscCallThrow(MatMultTranspose(cons->UnconMats[step], CForceVec2[step], fCon[step]));  //  F_con = U^T*f-U^T*K*C*C_v
    PetscCallThrow(VecScale(fCon[step], -1.0)); //Solving K*dx = -fint

    dofVersionUsed[step] = cons->dofVersion[step];
}

/// @brief Apply constraints to the force and stiffness matrix for the current staggered step
/// @param step staggered step number
void Physics::Constraint(size_t step){

    // stiffness matrix
    if (dofVersionUsed[step] != cons->dofVersion[step]){
        PetscCallThrow(MatProductCreate(K[step],cons->UnconMats[step],NULL, &KCon[step]));
        PetscCallThrow(MatProductSetType(KCon[step], MATPRODUCT_PtAP));
        PetscCallThrow(MatProductSetAlgorithm(KCon[step], MATPRODUCTALGORITHMDEFAULT));
        PetscCallThrow(MatProductSetFill(KCon[step], 1.0));
        PetscCallThrow(MatProductSetFromOptions(KCon[step]));
        PetscCallThrow(MatProductSymbolic(KCon[step]));
    }
    PetscCallThrow(MatProductNumeric(KCon[step]));  // K_con = U^T*K*U

    //force vector
    if (dofVersionUsed[step] != cons->dofVersion[step]){
        PetscCallThrow(VecCreateMPI(PETSC_COMM_WORLD, cons->ConstrainedSizeLocal[step], PETSC_DECIDE, &fCon[step]));
    }
    PetscCallThrow(MatMult(cons->ConMats[step], cons->ConValVec[step], CForceVec[step]));   //  F_c = C*C_v
    PetscCallThrow(VecScale(CForceVec[step], -1.0));
    PetscCallThrow(MatMultAdd(K[step], CForceVec[step], f[step], CForceVec2[step]));        //  F_c'= f-K*C*C_v
    PetscCallThrow(MatMultTranspose(cons->UnconMats[step], CForceVec2[step], fCon[step]));  //  F_con = U^T*f-U^T*K*C*C_v
    PetscCallThrow(VecScale(fCon[step], -1.0)); //Solving K*dx = -fint

    dofVersionUsed[step] = cons->dofVersion[step];

    //debugging, save to file
    if (exportMatrices){
        PetscViewer MatViewer;
        std::string FName = "matC"+std::to_string(step)+".bin";
        PetscViewerBinaryOpen(PETSC_COMM_WORLD, FName.c_str(), FILE_MODE_WRITE, &MatViewer);
        PetscCallThrow(MatView(KCon[step], MatViewer));
        PetscCallThrow(PetscViewerDestroy(&MatViewer));
    }
}

/// @brief Apply boundary conditions on top of the current and old state vectors
void Physics::BoundaryConsToOld(){
    //u_old = U u + C*C_v;
    for (size_t i = 0; i < dofspace->maxSteps; i++){
        for (const auto& [dof, value] : cons->conDofToVal[i]){
            StateVectorsOld[i].Set(dof, value, INSERT_VALUES);
            StateVectors[i].Set(dof, value, INSERT_VALUES);
        }
    }
    for (size_t i = 0; i < dofspace->maxSteps; i++){
        StateVectors[i].AssemblyStart();
        StateVectorsOld[i].AssemblyStart();
    }
    for (size_t i = 0; i < dofspace->maxSteps; i++){
        StateVectors[i].AssemblyEnd();
        StateVectorsOld[i].AssemblyEnd();
    }
    for (size_t i = 0; i < dofspace->maxSteps; i++){
        StateVectors[i].SyncStart(INSERT_VALUES);
        StateVectorsOld[i].SyncStart(INSERT_VALUES);
    }
    for (size_t i = 0; i < dofspace->maxSteps; i++){
        StateVectors[i].SyncEnd(INSERT_VALUES);
        StateVectorsOld[i].SyncEnd(INSERT_VALUES);
    }
    
    
}

void Physics::ResetState(){
    //reset state vectors to old values
    for (size_t i = 0; i < dofspace->maxSteps; i++){
        StateVectors[i].Set(StateVectorsOld[i]);
        dStateVectors[i].Set(dStateVectorsOld[i]);
        ddStateVectors[i].Set(ddStateVectorsOld[i]);
    }
    for (size_t i = 0; i < nModels; i++){
        Models[i]->ResetStep();
    }
}

/// @brief Adds iterative increment dx to the state vector, and updates time derivatives
/// @param step Staggered step
/// @param dx State vector increment
void Physics::UpdateState(size_t step, Vec& dx){
    // u = u_old + U dx + C*C_v;
    PetscCallThrow(MatMultAdd(cons->UnconMats[step], dx, StateVectors[step].DataVector, StateVectors[step].DataVector));
    PetscCallThrow(MatMultAdd(cons->ConMats[step], cons->ConValVec[step], StateVectors[step].DataVector, StateVectors[step].DataVector));
    
    StateVectors[step].AssemblyStart();
    StateVectors[step].AssemblyEnd();

    StateVectors[step].SyncStart(INSERT_VALUES);
    StateVectors[step].SyncEnd(INSERT_VALUES);
    UpdateVelAcc();
}

/// @brief Provides the requested nodal data for exporting results
/// @param NodalQuantityToSave input: Quantity to provide, identified by its string
/// @param EGroupIndex input: Index of the element group for which this quantity should be provided
/// @param NodeData output: values of quantity (should be pre-sized correctly)
void Physics::GetNodalDataToSave(std::string NodalQuantityToSave, size_t EGroupIndex, std::vector<std::vector<double>>& NodeData){
    bool DataSaved = false;

    size_t dofInd, dofStep;
    if (dofspace->hasDofType(NodalQuantityToSave, dofInd, dofStep)){ //directly get from dofspace
        size_t nNodes = mesh->ElementGroups[EGroupIndex].NNodes_per_elem;
        std::vector<size_t> NodeIdx(nNodes);
        std::vector<Eigen::RowVectorXd> NExport(mesh->ElementGroups[EGroupIndex].BaseElem->NExport.size()); 
        for (size_t i = 0; i < NExport.size(); i++) NExport[i].resize(nNodes);
        std::vector<size_t> DOFIdx(nNodes);
        Eigen::VectorXd NodeValues(nNodes); 
        for (size_t i = 0; i < mesh->ElementGroups[EGroupIndex].NElems; i++){
            mesh->getExportShape(EGroupIndex, i, NodeIdx, NExport);

            mesh->GetNodesForElem(NodeIdx, EGroupIndex, i);
            dofspace->getDofForNodes(NodeIdx, dofInd, DOFIdx);
            StateVectors[dofStep].GetValues(DOFIdx, NodeValues);
            for (size_t j = 0; j < NExport.size(); j++){
                NodeData[i][j] = NExport[j]*NodeValues;
            }
        }
        DataSaved = true;
    } else {
        size_t iModel = 0;
        while (DataSaved == false && iModel<nModels){
            DataSaved = Models[iModel]->SaveData("Nodes", NodalQuantityToSave, EGroupIndex, NodeData);
            iModel += 1;
        }
    }
    if (DataSaved == false){
        throw std::invalid_argument("Data of type \"" + NodalQuantityToSave + "\" is not saveable\n");
    }    
}

/// @brief Provides the requested integration-point-level data for exporting results
/// @param IPQuantityToSave input: Identifier for the integration-point quantity to save
/// @param EGroupIndex input: Index of the element group for which this quantity should be provided
/// @param IPData output: values of quantity (should be pre-sized correctly)
void Physics::GetIPDataToSave(std::string IPQuantityToSave, size_t EGroupIndex, std::vector<std::vector<double>>& IPData){
    size_t iModel = 0;
    bool DataSaved = false;
    while (DataSaved == false && iModel<nModels){
            DataSaved = Models[iModel]->SaveData("ip", IPQuantityToSave, EGroupIndex, IPData);
            iModel += 1;
        }
    if (DataSaved == false){
        throw std::invalid_argument("Data of type \"" + IPQuantityToSave + "\" is not saveable\n");
    }  
}

/// @brief updates and syncs between cores the time derivatives
void Physics::UpdateVelAcc(){
    for (size_t i = 0; i < dofspace->maxSteps; i++){
        timeScheme->UpdateVelAcc(StateVectors[i], dStateVectors[i], ddStateVectors[i], StateVectorsOld[i], dStateVectorsOld[i], ddStateVectorsOld[i]);
        dStateVectors[i].AssemblyStart();
        dStateVectors[i].AssemblyEnd();
        ddStateVectors[i].AssemblyStart();
        ddStateVectors[i].AssemblyEnd();
        dStateVectors[i].SyncStart(INSERT_VALUES);
        dStateVectors[i].SyncEnd(INSERT_VALUES);
        ddStateVectors[i].SyncStart(INSERT_VALUES);
        ddStateVectors[i].SyncEnd(INSERT_VALUES);
    }
}

/// @brief Commit history-dependent data
/// @param CommitType // TIMEDEP_COMMIT_TYPE = 1   or   PATHDEP_COMMIT_TYPE = 2;
void Physics::Commit(int CommitType){
    for (size_t i = 0; i < nModels; i++){
        Models[i]->Commit(CommitType);
    }

    if (CommitType == CommitTypes::TIMEDEP_COMMIT_TYPE){
        UpdateVelAcc();
        for (size_t i = 0; i < dofspace->maxSteps; i++){
            StateVectorsOld[i].Set(StateVectors[i]);
            dStateVectorsOld[i].Set(dStateVectors[i]);
            ddStateVectorsOld[i].Set(ddStateVectors[i]);

            StateVectorsOld[i].AssemblyStart();
            StateVectorsOld[i].AssemblyEnd();
            StateVectorsOld[i].SyncStart(INSERT_VALUES);
            StateVectorsOld[i].SyncEnd(INSERT_VALUES);

            dStateVectorsOld[i].AssemblyStart();
            dStateVectorsOld[i].AssemblyEnd();
            dStateVectorsOld[i].SyncStart(INSERT_VALUES);
            dStateVectorsOld[i].SyncEnd(INSERT_VALUES);

            ddStateVectorsOld[i].AssemblyStart();
            ddStateVectorsOld[i].AssemblyEnd();
            ddStateVectorsOld[i].SyncStart(INSERT_VALUES);
            ddStateVectorsOld[i].SyncEnd(INSERT_VALUES);
        }

        Append_Timeseries();

        #ifdef LOCALENVIRONMENT
            Vis_Plots->Plot("step", 0);
        #endif
    }
}

void Physics::Init_Timeseries(){
    TimeDataTypes.resize(0);

    size_t nData = 3;
    TimeDataTypes.push_back("time");
    TimeDataTypes.push_back("dt");
    TimeDataTypes.push_back("tstep");

    for (size_t i = 0; i < nModels; i++){
        std::vector<std::string> dataNames(0);
        dataNames.resize(0);
        nData += Models[i]->hasTimeData(dataNames);
        for (size_t j = 0; j < dataNames.size(); j++){
            TimeDataTypes.push_back(dataNames[j]);
        }
    }

    TimeData.resize(nData);
    for (size_t i = 0; i < nData; i++){
        TimeData[i].resize(0);
        TimeData[i].push_back(0.0);
    }
}

double Physics::GetTimeData(std::string DataName){
    double tData = 0;
    size_t nData = 0;
    for (size_t i = 0; i < nModels; i++){
        std::vector<std::string> dataNames(0);
        nData = Models[i]->hasTimeData(dataNames);
        for (size_t j = 0; j < nData; j++){
            if (DataName == dataNames[j]){
                std::vector<double> dataVals;
                Models[i]->GetTimeData(dataVals);
                tData = dataVals[j];
                return tData;
            }
        }
    }
    return 0;
}

void Physics::Append_Timeseries(){
    size_t iData = 0;
    TimeData[iData].push_back(time+timeScheme->dt);

    iData = 1;
    TimeData[iData].push_back(timeScheme->dt);

    iData = 2;
    TimeData[iData].push_back(TimeData[iData].size());

    for (size_t i = 0; i < nModels; i++){
        std::vector<double> dataVals;

        dataVals.resize(0);
        Models[i]->GetTimeData(dataVals);
        for (size_t j = 0; j < dataVals.size(); j++){
            iData += 1;
            TimeData[iData].push_back(dataVals[j]);
        }
    }
}
