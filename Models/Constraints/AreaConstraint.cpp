#include "AreaConstraint.h"
#include "../../Physics/physics.h"

void Register_AreaConstraintsModel() {
    ModelNames.push_back("AreaConstraints");
    ModelCreators.push_back(New_AreaConstraintsModel);
}

BaseModel* New_AreaConstraintsModel(Physics& My_Physics, std::string MyNameIn) {
    return new AreaConstraintsModel(My_Physics, MyNameIn);
}

AreaConstraintsModel::AreaConstraintsModel(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "AreaConstraints";
};

AreaConstraintsModel::~AreaConstraintsModel() {

}

void AreaConstraintsModel::load(inputData& inputs, SaveDataFile& data) {
    Setup(inputs);
}

void AreaConstraintsModel::save(SaveDataFile& data) {
    // Save implementation
}

void AreaConstraintsModel::init(inputData& inputs) {
    Setup(inputs);
    std::vector<size_t> UniqueNodes = mesh->ElementGroups[ElementGroupIndex].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes, dofTypes);
}

void AreaConstraintsModel::Setup(inputData& inputs) {
    std::string NodeGroupName; inputs.GetRequired(NodeGroupName, {"Models", MyName.c_str(), "NodeGroup"});
    ElementGroupIndex = physics->mesh->GetNodeGroupIdx(NodeGroupName);
    
    // Get dofs
    inputs.GetRequired(dofNames, {"Models",MyName,"dofs"});
    inputs.GetRequired(conVals, {"Models",MyName,"values"});

    dofs->getDofTypesSteps(dofNames, dofTypes, dofSteps);

    xLims.resize(2); yLims.resize(2); zLims.resize(2);
    xLims[0] = -1e99; xLims[1] = 1e99;
    inputs.GetOptional(xLims, {"Models", MyName, "xLims"});

    yLims[0] = -1e99; yLims[1] = 1e99;
    inputs.GetOptional(yLims, {"Models", MyName, "yLims"});

    zLims[0] = -1e99; zLims[1] = 1e99;
    inputs.GetOptional(zLims, {"Models", MyName, "zLims"});

    kDummy = 1.0e16;
    inputs.GetOptional(kDummy, {"Models", MyName, "kDummy"});
}

void AreaConstraintsModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step) {
    for (size_t conDOF = 0; conDOF < dofTypes.size(); conDOF++){
        if (step == dofSteps[conDOF]){
            std::vector<size_t> UniqueNodes = mesh->NodeGroups[ElementGroupIndex].Nodes;
            size_t nNodes = UniqueNodes.size();
            std::vector<size_t> NodeDofs(nNodes);
            dofs->getDofForNodes(UniqueNodes, dofTypes[conDOF], NodeDofs);
            std::vector<double> X(nNodes), Y(nNodes), Z(nNodes); 
            if (dim==2){
                mesh->GetCoordsForNodes(X, Y, UniqueNodes);
            } else {
                mesh->GetCoordsForNodes(X, Y, Z, UniqueNodes);
            }

            for (size_t i = 0; i < nNodes; i++){
                double x = X[i], y = Y[i], z = 0.0; 
                if (dim==3){z = Z[i];};
                if (x > xLims[0] && x < xLims[1] && y > yLims[0] && y < yLims[1] && z > zLims[0] && z < zLims[1]){
                    cons->AddConstraint(dofSteps[conDOF], NodeDofs[i], conVals[conDOF]);
               }
            }
        }
    }
}

