#include "BaseLinSolver.h"

BaseLinSolver::BaseLinSolver(){

}

BaseLinSolver::~BaseLinSolver(){

}

/// @brief Initializes linear solver based on input data file
/// @param inputs Input data file
/// @param K matrix to use during solving
void BaseLinSolver::init(inputData& inputs, Mat& K){

}

/// @brief Solve linear system of equations
/// @param K input: system matrix
/// @param f input: force vector
/// @param dx output: state increment
void BaseLinSolver::solve(Mat& K, Vec& f, Vec& dx){

}