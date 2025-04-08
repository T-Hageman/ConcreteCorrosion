#ifndef GMRESLINSOLVER_H
#define GMRESLINSOLVER_H

#include <vector>

#include "BaseLinSolver.h"

/// @brief GMRES direct linear solver
class GMRESLinSolver: public BaseLinSolver {
    public:
        std::string MyName;

        GMRESLinSolver();
        ~GMRESLinSolver();

        void init(inputData& inputs, Mat& K);
        void solve(Mat& K, Vec& f, Vec& dx);
    protected:

    private:
        KSP ksp;    //PETSC solver object
        PC pc;  //actual linear solver
        Mat PardisoMat;
        bool firstSolve;

        std::vector<PetscInt> IOptionIdx, IOptionVal;   //combination of input option index and value to set it to
};
#endif