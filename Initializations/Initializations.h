#ifndef INITIALISATIONS_H
#define INITIALISATIONS_H

#include "../InputsOutputs/inputData.h"
#include "../Physics/physics.h"

void Initialize(inputData& inputs, Physics& MyPhysics);
void Initialize_Zero(Physics& MyPhysics);
void Initialize_ElectroChemistry(inputData& inputs, Physics& MyPhysics);


#endif