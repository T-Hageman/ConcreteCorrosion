#ifndef PARDISOLINSOLVER_H
#define PARDISOLINSOLVER_H

#include <vector>

#include "BaseLinSolver.h"

/// @brief PARDISO direct linear solver. For input parameters see 
/// https://pardiso-project.org/manual/manual.pdf or https://www.intel.com/content/www/us/en/develop/documentation/onemkl-developer-reference-c/top/sparse-solver-routines/onemkl-pardiso-parallel-direct-sparse-solver-iface/pardiso-iparm-parameter.html"
class PARDISOLinSolver: public BaseLinSolver {
    public:
        std::string MyName;

        PARDISOLinSolver();
        ~PARDISOLinSolver();

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