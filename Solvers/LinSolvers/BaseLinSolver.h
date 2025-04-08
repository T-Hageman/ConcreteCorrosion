#ifndef LINSOLVER_H
#define LINSOLVER_H

#include "../../InputsOutputs/inputData.h"

#include <petsc.h>
#include <iostream>

/// @brief Linear solver to solve K*dx = -f
class BaseLinSolver {
    public:
        std::string MyName;

        BaseLinSolver();
        virtual ~BaseLinSolver();

        virtual void init(inputData& inputs, Mat& K);
        virtual void solve(Mat& K, Vec& f, Vec& dx);

        bool Update_K = true;
    protected:

    private:




};
#endif