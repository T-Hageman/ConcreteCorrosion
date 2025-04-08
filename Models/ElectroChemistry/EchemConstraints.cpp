#include "EchemConstraints.h"
#include "../../Physics/physics.h"

void Register_EchemConstraints(){
    ModelNames.push_back("Electrochemistry/EchemConstraints");
    ModelCreators.push_back(New_EchemConstraints);
}
BaseModel* New_EchemConstraints(Physics& My_Physics, std::string MyNameIn){
    return new EchemConstraints(My_Physics, MyNameIn);
}

EchemConstraints::EchemConstraints(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "Electrochemistry/EchemConstraints";
}

EchemConstraints::~EchemConstraints(){

}

void EchemConstraints::init(inputData& inputs){
    Setup(inputs);


    dofs->AddDofs(mesh->NodeGroups[NodeGroupIndex_C].Nodes, dofTypes_C);
    dofs->AddDofs(mesh->NodeGroups[NodeGroupIndex_E].Nodes, dofType_ePot);
}

void EchemConstraints::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
}

void EchemConstraints::save(SaveDataFile& data){

}

void EchemConstraints::Commit(int CommitType){

}

void EchemConstraints::Setup(inputData& inputs){
    std::string GroupName_E, GroupName_C;
    inputs.GetRequired(GroupName_C, {"Models",MyName,"NodeGroup_C"});
    inputs.GetRequired(GroupName_E, {"Models",MyName,"NodeGroup_E"});

    NodeGroupIndex_C = mesh->GetNodeGroupIdx(GroupName_C);
    NodeGroupIndex_E = mesh->GetNodeGroupIdx(GroupName_E);
    
    // Get dofs
    inputs.GetRequired(Species, {"Models",MyName,"Species"});
    nSpecies = Species.size();

    ChemPotBased = false;
    inputs.GetOptional(ChemPotBased, {"properties","Species","UseChemPotentials"});

    std::vector<size_t> dofSteps;
    dofs->getDofTypesSteps(Species, dofTypes_C, dofSteps);
    Step_C = dofSteps[0];
    for (size_t i = 1; i < nSpecies; i++){
        if (Step_C!=dofSteps[i]){
            throw std::invalid_argument(ModelName+" requires "+ Species[0] + " and " + Species[i] + " to be in the same solver step\n");
        }
    }
    size_t EStep;
    dofs->getDofTypesSteps(DofNames_ePot, dofType_ePot, EStep);
    if (EStep != Step_C){
        throw std::invalid_argument(ModelName+" requires "+ Species[0] + " and " + DofNames_ePot + " to be in the same solver step\n");
    }
    inputs.GetRequired(E0, {"Models",MyName.c_str(),"ePot"});

    C0.resize(nSpecies);
    z.resize(nSpecies);

    int i_from_neutrality = -1;
    for (size_t i = 0; i < nSpecies; i++){
        inputs.GetRequired(z[i], {"properties","Species",Species[i],"z"});
        std::string C0Type = inputs.GetType({"properties","Species",Species[i].c_str(),"C0"});

        if (C0Type == "Number"){
            inputs.GetRequired(C0[i], {"properties","Species",Species[i].c_str(),"C0"});
            if (ChemPotBased && C0[i]<std::exp(-30.0)){
                C0[i] = std::exp(-30.0);
            }
        } else if (C0Type == "String") {
            if (i_from_neutrality<0){
                i_from_neutrality = i;
            } else {
                throw std::invalid_argument("Multiple species are designated as initialized from electroneutrality\n");
            }
        } else {
            inputs.GetRequired(C0[i], {"properties","Species",Species[i].c_str(),"C0"}); //generate error here
        }
    }
    if (i_from_neutrality>=0){
        double Z_Filler = 0;
        for (size_t i = 0; i < nSpecies; i++){
            int ii = i;
            if (ii!=i_from_neutrality){
                Z_Filler += z[i]*C0[i];
            }
        }
        C0[i_from_neutrality] = -Z_Filler/z[i_from_neutrality];
        if (C0[i_from_neutrality]<1.0e-20){
            throw std::invalid_argument("Initial concentrations do not follow from electroneutrality\n");
        }
    }
}

void EchemConstraints::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    //double t = physics->time;
    if (step == Step_C){
        for (size_t i = 0; i < nSpecies; i++){
            std::vector<size_t> dofIDs(mesh->NodeGroups[NodeGroupIndex_C].NNodes);
            dofs->getDofForNodes(mesh->NodeGroups[NodeGroupIndex_C].Nodes, dofTypes_C[i], dofIDs);
            double cCon;
            if (ChemPotBased){
                cCon = std::log(C0[i]);
            } else {
                cCon = C0[i];
            }
            //std::cout << C0[i] << "->" << cCon << "\n";
            cons->AddConstraint(Step_C, dofIDs, cCon);
        }
        std::vector<size_t> dofIDs(mesh->NodeGroups[NodeGroupIndex_E].NNodes);
        dofs->getDofForNodes(mesh->NodeGroups[NodeGroupIndex_E].Nodes, dofType_ePot, dofIDs);
        cons->AddConstraint(Step_C, dofIDs, E0);  
    }
}