#include "Constrainer.h"

Constrainer::Constrainer(DofSpace* dofspace_in){
    PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));

    dofspace = dofspace_in;

    conDofToVal.resize(dofspace->maxSteps);
    ConMats.resize(dofspace->maxSteps);
    UnconMats.resize(dofspace->maxSteps);
    ConValVec.resize(dofspace->maxSteps);
    dofsChanged.resize(dofspace->maxSteps);
    ConStartNumber.resize(dofspace->maxSteps);
    dofVersion.resize(dofspace->maxSteps);
    ConstrainedSizeLocal.resize(dofspace->maxSteps);

    matInitSize.resize(dofspace->maxSteps);
    for (size_t i = 0; i < dofspace->maxSteps; i++){
        matInitSize[i] = 0;
        PetscCallThrow(MatCreateAIJ(PETSC_COMM_WORLD, 0, 0, 0, 0, 1, NULL, 1, NULL, &ConMats[i]));
        PetscCallThrow(MatCreateAIJ(PETSC_COMM_WORLD, 0, 0, 0, 0, 1, NULL, 1, NULL, &UnconMats[i]));

        dofsChanged[i] = true;
        dofVersion[i] = 0;
    }
    
};

Constrainer::~Constrainer(){

};

/// @brief Empties the constraints
/// @param step staggered step
void Constrainer::SetZero(size_t step){
    //size_t reserves = conDofToVal.size();
    conDofToVal[step].clear();
}

/// @brief Adds constrained dofs (and checks for duplicates)
/// @param step Staggered step for which to constrain
/// @param dof Degree of freedom being constrained
/// @param value Value being applied as constraint
void Constrainer::AddConstraint(size_t step, size_t dof, double value){
    auto it = conDofToVal[step].find(dof);
    if (it != conDofToVal[step].end()){
        if (it->second != value){
            throw std::invalid_argument("Degree of freedom "+ std::to_string(dof) + " is doubly constrained with different values\n");
        }
    } else {
        conDofToVal[step].insert({dof, value});
    }
}

/// @brief Adds constrained dofs (and checks for duplicates)
/// @param step Staggered step for which to constrain
/// @param dof Degree of freedom being constrained
/// @param value Value being applied as constraint
void Constrainer::AddConstraint(size_t step, std::vector<size_t> &dof, double value){
    for (size_t i = 0; i < dof.size(); i++){
        AddConstraint(step, dof[i], value);
    }
}

/// @brief Assembles constrain matrices (does not apply these constraints, just assembles structures required for it)
/// @param step Staggered steps to constrain for
/// @param State State vector for current step
void Constrainer::Assemble(size_t step, GhostVector& State){
    size_t numCons = conDofToVal[step].size();
    size_t numDofs = dofspace->LocalDofNumRange[step][1]-dofspace->LocalDofNumRange[step][0];
    size_t numUncon= numDofs-numCons;

    if (dofsChanged[step]){
        NumConsGlobal.resize(size);
        for (int i = 0; i < size; i++){
            NumConsGlobal[i] = 0; 
        }
        NumConsGlobal[rank] = numCons;
        MPI_Allreduce(MPI_IN_PLACE, NumConsGlobal.data(), size, my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
        size_t StartNumCon = 0; 
        for (int i = 0; i < rank; i++){
            StartNumCon += NumConsGlobal[i];
        }
        ConStartNumber[step] = StartNumCon;
        ConstrainedSizeLocal[step] = numUncon;
        size_t StartNumUncon = dofspace->LocalDofNumRange[step][0]-StartNumCon;

        if (matInitSize[step] == numUncon){
            PetscCallThrow(MatZeroEntries(ConMats[step]));
            PetscCallThrow(MatZeroEntries(UnconMats[step]));
        } else {
            PetscCallThrow(MatDestroy(&UnconMats[step]));
            PetscCallThrow(MatDestroy(&ConMats[step]));
            PetscCallThrow(MatCreateAIJ(PETSC_COMM_WORLD, numDofs, numUncon,  PETSC_DECIDE, PETSC_DECIDE, 1, NULL, 1, NULL, &UnconMats[step]));
            PetscCallThrow(MatCreateAIJ(PETSC_COMM_WORLD, numDofs, numCons,  PETSC_DECIDE, PETSC_DECIDE, 1, NULL, 1, NULL, &ConMats[step]));
            matInitSize[step] = numUncon;
        }

        for (size_t dof = dofspace->LocalDofNumRange[step][0]; dof < dofspace->LocalDofNumRange[step][1]; dof++){
            bool hasCon = conDofToVal[step].count(dof);
            if (hasCon) {
                PetscCallThrow(MatSetValue(ConMats[step], dof, StartNumCon, 1, INSERT_VALUES));
                StartNumCon += 1;
            } else {
                PetscCallThrow(MatSetValue(UnconMats[step], dof, StartNumUncon, 1, INSERT_VALUES));
                StartNumUncon += 1;
            }
        }
        VecCreateMPI(PETSC_COMM_WORLD, numCons, PETSC_DECIDE, &ConValVec[step]);
        dofsChanged[step] = false;
        dofVersion[step] += 1;
    }

    VecZeroEntries(ConValVec[step]);
    size_t iLoc = ConStartNumber[step];
    for (const auto& [dof, value] : conDofToVal[step]){
        double conVal_dx = value-State.GetValue(dof);
        VecSetValue(ConValVec[step], iLoc, conVal_dx, INSERT_VALUES);
        iLoc += 1;
    }
    
    PetscCallThrow(MatAssemblyBegin(ConMats[step], MAT_FINAL_ASSEMBLY));
    PetscCallThrow(MatAssemblyBegin(UnconMats[step], MAT_FINAL_ASSEMBLY));
    PetscCallThrow(VecAssemblyBegin(ConValVec[step]));
}

/// @brief End of constraint object assembly
/// @param step 
void Constrainer::AssembleEnd(size_t step){
    PetscCallThrow(MatAssemblyEnd(ConMats[step], MAT_FINAL_ASSEMBLY));
    PetscCallThrow(MatAssemblyEnd(UnconMats[step], MAT_FINAL_ASSEMBLY));
    PetscCallThrow(VecAssemblyEnd(ConValVec[step]));
}