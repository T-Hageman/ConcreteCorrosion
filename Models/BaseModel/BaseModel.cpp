#include "BaseModel.h"
#include "../../Physics/physics.h"

void Register_BaseModel(){
    ModelNames.push_back("BaseModel");
    ModelCreators.push_back(New_BaseModel);
}

BaseModel* New_BaseModel(Physics& My_Physics, std::string MyNameIn){
    return new BaseModel(My_Physics, MyNameIn);
}

BaseModel::BaseModel(Physics& My_Physics, std::string MyNameIn){
    ModelName = "BaseModel";
    MyName = MyNameIn;

    physics = &My_Physics;
    mesh = physics->mesh;
    dofs = physics->dofspace;

    dim = mesh->dim;

    MPI_Comm_size(PETSC_COMM_WORLD,&MPI_size);
    MPI_Comm_rank(PETSC_COMM_WORLD,&MPI_rank);
};

BaseModel::~BaseModel(){

};

void BaseModel::ResetStep(){
    
}

void BaseModel::load(inputData& inputs, SaveDataFile& data){

}

void BaseModel::save(SaveDataFile& data){

}

void BaseModel::init(inputData& inputs){
    
}

void BaseModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    
}

bool BaseModel::SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data){
    return false;
}

void BaseModel::Commit(int CommitType){

}

size_t BaseModel::hasTimeData(std::vector<std::string>& DataNames){
    return 0;
}

void BaseModel::GetTimeData(std::vector<double>& DataValues){

}