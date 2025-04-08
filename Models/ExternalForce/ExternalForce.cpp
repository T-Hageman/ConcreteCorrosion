#include "ExternalForce.h"
#include "../../Physics/physics.h"

void Register_ExternalForceModel(){
    ModelNames.push_back("ExternalForce");
    ModelCreators.push_back(New_ExternalForceModel);
}

BaseModel* New_ExternalForceModel(Physics& My_Physics, std::string MyNameIn){
    return new ExternalForceModel(My_Physics, MyNameIn);
}

ExternalForceModel::ExternalForceModel(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "ExternalForce";
};

ExternalForceModel::~ExternalForceModel(){

};

void ExternalForceModel::init(inputData& inputs){
    Setup(inputs);

    std::vector<size_t> UniqueNodes = mesh->ElementGroups[ElementGroupIndex].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes, dofTypes);
}

void ExternalForceModel::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
}

void ExternalForceModel::Setup(inputData& inputs){
    std::string EGroupName; inputs.GetRequired(EGroupName,{"Models",MyName,"ElementGroup"});
    ElementGroupIndex = mesh->GetElementGroupIdx(EGroupName);
    
    // Get and Add dofs
    inputs.GetRequired(dofNames,{"Models",MyName,"dofs"});
    inputs.GetRequired(ForceVals,{"Models",MyName,"Force"});
    dofs->getDofTypesSteps(dofNames, dofTypes, dofSteps);
}

void ExternalForceModel::save(SaveDataFile& data){

}

void ExternalForceModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    for (size_t i = 0; i < dofNames.size(); i++){
        if (dofSteps[i] == step){
            size_t ipcount = mesh->ElementGroups[ElementGroupIndex].BaseElem->ipcount;
            size_t nNodes = mesh->ElementGroups[ElementGroupIndex].NNodes_per_elem;
            Eigen::VectorXd F_el(nNodes);
            Eigen::MatrixXd coordsNodes(nNodes, 2);
            std::vector<Eigen::RowVectorXd> N(ipcount); for (size_t ip=0; ip<ipcount;  ip++) N[ip].resize(nNodes); 
            std::vector<Eigen::MatrixXd> G(ipcount); for (size_t ip=0; ip<ipcount;  ip++) G[ip].resize(1, nNodes); 
            std::vector<double> w(ipcount);
            std::vector<size_t> X(nNodes);
            std::vector<PetscInt> dofsX(nNodes);

            for (size_t el = 0; el < mesh->ElementGroups[ElementGroupIndex].NElems; el++){ 
                mesh->getShapeGrads(ElementGroupIndex, el, X, N, G, w);

                dofs->getDofForNodes(X, dofTypes[i], dofsX);
                F_el.setZero(); 
                for (size_t ip = 0; ip < ipcount; ip++){
                    F_el += -w[ip] * N[ip].transpose()*ForceVals[i];
                }
                VecAdd(f, dofsX, F_el);
            }
        }
    }
}