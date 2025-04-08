#include "ModelRegister.h"
#include <iostream>

std::vector<std::string> ModelNames; //vector of available model names
std::vector<std::function<BaseModel*(Physics&, std::string)>> ModelCreators; //pointer to model creation functions


/// @brief Registers the available physics models
void RegisterModels(){
    //Basic models
    Register_BaseModel();
    Register_ConstraintsModel();
    Register_ExternalForceModel();
    Register_WeakExternalModel();
    Register_WeakAreaModel();
    Register_AreaConstraintsModel();

    //electrochemistry
    Register_NernstPlanck();
    Register_EchemConstraints();
    Register_ElectroSurface();
}


/// @brief Creates the physics models
/// @param physics reference to physics object
/// @param ModelName characteristic model to be created
/// @param MyName unique identifier to name model
/// @return pointer to model (remember to destroy when not needed)
BaseModel* CreateModel(Physics& physics, std::string ModelName, std::string MyName){
    BaseModel* Model;

    size_t it = 0;
    while (true){
        if (ModelName == ModelNames[it]){
            Model = ModelCreators[it](physics, MyName);
            return Model;
        }
        it+= 1;
        if (it == ModelNames.size()){
            std::string validTypes = "";
            for (size_t i = 0; i < ModelNames.size(); i++){
                validTypes.append(ModelNames[i]);
                validTypes.append(", ");
            }
            throw std::invalid_argument("Model type "+ModelName+" not defined, valid options are: " + validTypes);
            return 0;
        }
    } 
}