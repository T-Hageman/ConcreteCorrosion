#ifndef MODELTYPEREGISTER_H
#define MODELTYPEREGISTER_H

#include "../InputsOutputs/inputData.h"
#include <string>

#include "BaseModel/BaseModel.h"
#include "Constraints/ConstraintsModel.h"
#include "ExternalForce/ExternalForce.h"
#include "ExternalForce/WeakExternal.h"
#include "ExternalForce/WeakArea.h"
#include "ElectroChemistry/NernstPlanck.h"
#include "ElectroChemistry/EchemConstraints.h"
#include "ElectroChemistry/ElectroSurface.h"
#include "Constraints/AreaConstraint.h"


extern std::vector<std::string> ModelNames;
extern std::vector<std::function<BaseModel*(Physics&, std::string)>> ModelCreators;

void RegisterModels();
BaseModel* CreateModel(Physics& physics, std::string ModelName, std::string MyName);

#endif