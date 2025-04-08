#include "Initializations.h"
#include <algorithm>

/// @brief Sets initial conditions based on input parameter "nonlinSolver.Initialization.Type"
/// @param inputs Input properties object
/// @param MyPhysics pointer to physics object
void Initialize(inputData& inputs, Physics& MyPhysics){
    std::string InitType; inputs.GetRequired(InitType,{"nonlinSolver","Initialization","Type"});

    if (InitType == "Zero") Initialize_Zero(MyPhysics);
    if (InitType == "ElectroChemistry") Initialize_ElectroChemistry(inputs, MyPhysics);


    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        MyPhysics.StateVectors[i].AssemblyStart();
        MyPhysics.dStateVectors[i].AssemblyStart();
        MyPhysics.ddStateVectors[i].AssemblyStart();
    }
    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        MyPhysics.StateVectors[i].AssemblyEnd();
        MyPhysics.dStateVectors[i].AssemblyEnd();
        MyPhysics.ddStateVectors[i].AssemblyEnd();
        MyPhysics.StateVectors[i].SyncStart(INSERT_VALUES);
        MyPhysics.dStateVectors[i].SyncStart(INSERT_VALUES);
        MyPhysics.ddStateVectors[i].SyncStart(INSERT_VALUES);
    } 
    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        MyPhysics.StateVectors[i].SyncEnd(INSERT_VALUES);
        MyPhysics.dStateVectors[i].SyncEnd(INSERT_VALUES);
        MyPhysics.ddStateVectors[i].SyncEnd(INSERT_VALUES);
    }
}

/// @brief Initializer to set properties to zero (called from Initialize, do not call directly)
/// @param MyPhysics pointer to physics object
void Initialize_Zero(Physics& MyPhysics){
    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        VecZeroEntries(MyPhysics.StateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.dStateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.ddStateVectors[i].DataVector);
    }
}

void Initialize_ElectroChemistry(inputData& inputs, Physics& MyPhysics){
    for (size_t i = 0; i < MyPhysics.dofspace->maxSteps; i++){
        VecZeroEntries(MyPhysics.StateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.dStateVectors[i].DataVector);
        VecZeroEntries(MyPhysics.ddStateVectors[i].DataVector);
    }


    bool ChemPotBased; inputs.GetRequired(ChemPotBased, {"properties","Species","UseChemPotentials"});
    size_t nSpecies;
    std::vector<std::string> Species; inputs.GetRequired(Species, {"properties","Species","Types"});
    std::vector<double> C0, z;
    int i_from_neutrality = -1;
    double Z_Filler = 0;


    nSpecies = Species.size();
    C0.resize(nSpecies);
    z.resize(nSpecies);
    for (size_t i = 0; i < nSpecies; i++){
        inputs.GetRequired(z[i], {"properties","Species",Species[i],"z"});
        std::string C0Type = inputs.GetType({"properties","Species",Species[i].c_str(),"C0"});
        if (C0Type == "Number"){
            inputs.GetRequired(C0[i], {"properties","Species",Species[i],"C0"});
            if (ChemPotBased && C0[i]<std::exp(-30.0)){
                C0[i] = std::exp(-30.0);
            }
            Z_Filler += z[i]*C0[i];
        } else {
            if (i_from_neutrality<0){
                i_from_neutrality = i;
            } else {
                throw std::invalid_argument("Multiple species are designated as initialized from electroneutrality\n");
            }
        }
    }
    if (i_from_neutrality>=0){
        C0[i_from_neutrality] = -Z_Filler/z[i_from_neutrality];
        if (C0[i_from_neutrality]<1.0e-20){
            throw std::invalid_argument("Initial concentrations do not follow from electroneutrality\n");
        }
    }

    std::string InteriorGroup; inputs.GetRequired(InteriorGroup,{"nonlinSolver","Initialization","InteriorGroup"});
    size_t InteriorGroupIdx = MyPhysics.mesh->GetNodeGroupIdx(InteriorGroup);

    std::vector<size_t> DofTypes(nSpecies), dofSteps(nSpecies);
    MyPhysics.dofspace->getDofTypesSteps(Species, DofTypes, dofSteps);
    std::vector<size_t> Nodes = MyPhysics.mesh->NodeGroups[InteriorGroupIdx].Nodes;
    std::vector<double> NodeValues(Nodes.size());
    std::vector<PetscInt> NodeDofs(Nodes.size());

    for (size_t i = 0; i < nSpecies; i++){
        MyPhysics.dofspace->getDofForNodes(Nodes, DofTypes[i], NodeDofs);
        double cCon;
        if (ChemPotBased){
            cCon = std::log(C0[i]);
        } else {
            cCon = C0[i];
        }
        
        for (size_t j = 0; j < Nodes.size(); j++){
            NodeValues[j] = cCon;
        }
        MyPhysics.StateVectors[dofSteps[i]].Set(NodeDofs, NodeValues, INSERT_VALUES);
    }
}
