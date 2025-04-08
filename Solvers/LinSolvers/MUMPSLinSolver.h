#ifndef MUMPSLINSOLVER_H
#define MUMPSLINSOLVER_H

#include <vector>

#include "BaseLinSolver.h"


/// @brief MUMPS direct linear solver. For input parameters see https://graal.ens-lyon.fr/MUMPS/doc/userguide_5.5.1.pdf
class MUMPSLinSolver: public BaseLinSolver {
    public:
        std::string MyName;

        MUMPSLinSolver();
        ~MUMPSLinSolver();

        void init(inputData& inputs, Mat& K);
        void solve(Mat& K, Vec& f, Vec& dx);
    protected:

    private:
        KSP ksp;    //PETSC solver object
        PC pc;  //actual linear solver

        std::vector<PetscInt> IOptionIdx, IOptionVal;   //combination of input option index and value to set it to

};
#endif