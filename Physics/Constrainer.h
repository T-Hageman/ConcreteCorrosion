#ifndef CONSTRAINER_H
#define CONSTRAINER_H

#include <vector>
#include <Eigen/Dense>
#include <petscvec.h>
#include <map>

#include "../utility/utility.h"
#include "DofSpace.h"

/// @brief Handles nodal constraints
class Constrainer {
    public:
        Constrainer(DofSpace* dofspace_in);
        ~Constrainer();

        void SetZero(size_t step);
        void AddConstraint(size_t step, size_t dof, double value);
        void AddConstraint(size_t step, std::vector<size_t> &dof, double value);
        void Assemble(size_t step, GhostVector& State);
        void AssembleEnd(size_t step);

        DofSpace* dofspace; //reference to dofspace object

        std::vector<std::map<size_t, double>> conDofToVal;  //saves constraints being applied to each dof
        std::vector<Mat> ConMats, UnconMats;    //reorder matrices for applying constraints
        std::vector<size_t> ConstrainedSizeLocal;   //size of constrained problem
        std::vector<Vec> ConValVec;     //values of constrained nodes
        std::vector<bool> dofsChanged;  //should the data structures be re-calculated
        std::vector<size_t> dofVersion; //version of dof matrices

    private:
        PetscInt rank, size;
        std::vector<size_t> matInitSize;
        std::vector<size_t> NumConsGlobal;
        std::vector<size_t> ConStartNumber;
};

#endif