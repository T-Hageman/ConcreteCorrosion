#ifndef LINSOLVETYPEREGISTER_H
#define LINSOLVETYPEREGISTER_H

#include "../../InputsOutputs/inputData.h"
#include <string>

#include "BaseLinSolver.h"

BaseLinSolver* CreateLinSolver(std::string SolverType);

#endif