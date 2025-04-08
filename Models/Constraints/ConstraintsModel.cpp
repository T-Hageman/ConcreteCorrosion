#include "ConstraintsModel.h"
#include "../../Physics/physics.h"

void Register_ConstraintsModel(){
    ModelNames.push_back("Constraints");
    ModelCreators.push_back(New_ConstraintsModel);
}

BaseModel* New_ConstraintsModel(Physics& My_Physics, std::string MyNameIn){
    return new ConstraintsModel(My_Physics, MyNameIn);
}

ConstraintsModel::ConstraintsModel(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "Constraints";
};

ConstraintsModel::~ConstraintsModel(){

};

void ConstraintsModel::init(inputData& inputs){
    Setup(inputs);
    dofs->AddDofs(mesh->NodeGroups[NodeGroupIndex].Nodes, dofTypes);
}

void ConstraintsModel::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);

}

void ConstraintsModel::save(SaveDataFile& data){
    
}

void ConstraintsModel::Setup(inputData& inputs){
    std::string NodeGroupName; inputs.GetRequired(NodeGroupName, {"Models", MyName.c_str(), "NodeGroup"});
    NodeGroupIndex = physics->mesh->GetNodeGroupIdx(NodeGroupName);
    
    // Get dofs
    inputs.GetRequired(dofNames, {"Models",MyName,"dofs"});
    inputs.GetRequired(conVals, {"Models",MyName,"values"});

    dofs->getDofTypesSteps(dofNames, dofTypes, dofSteps);
}

void ConstraintsModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    for (size_t i = 0; i < dofNames.size(); i++){
        if (dofSteps[i] == step){
            std::vector<size_t> dofIDs(mesh->NodeGroups[NodeGroupIndex].NNodes);
            dofs->getDofForNodes(mesh->NodeGroups[NodeGroupIndex].Nodes, dofTypes[i], dofIDs);
            cons->AddConstraint(step, dofIDs, conVals[i]);
        }
    }
}