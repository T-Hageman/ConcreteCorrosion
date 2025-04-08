#ifndef DOFSPACE_H
#define DOFSPACE_H

#include <vector>
#include <Eigen/Dense>
#include <petscvec.h>

#include "../utility/utility.h"
#include "../mesh/mesh.h"
#include "../InputsOutputs/inputData.h"

/// @brief Class which handles degree of freedom numbering across all staggered steps
class DofSpace {
    public:
        DofSpace(inputData& inputs);
        ~DofSpace();

        void save(SaveDataFile& data);
        void load(SaveDataFile& data);

        void getDofTypesSteps(std::vector<std::string> names, std::vector<size_t> &DofTypes_out, std::vector<size_t> &DofSteps_out);
        void getDofTypesSteps(std::string name, size_t &DofType_out, size_t &DofStep_out);
        bool hasDofType(std::string name, size_t &DofType_out, size_t &DofStep_out);
        bool hasDofType(std::string name);
        void AddDofs(std::vector<size_t> &Nodes, std::vector<size_t> &DofInd);
        void AddDofs(std::vector<size_t> &Nodes, size_t DofInd);
        void AddDofs(size_t Node, size_t DofInd);
        void SyncDofs(Mesh* mesh);
        void printStats(size_t curStep);
        void ExportDofSpace();

        void getDofForNodes(size_t &Node, size_t DofInd, PetscInt &Dof_out);
        void getDofForNodes(std::vector<size_t> &Nodes, size_t DofInd, std::vector<size_t> &Dofs_out);
        void getDofForNodes(std::vector<size_t> &Nodes, size_t DofInd, std::vector<PetscInt> &Dofs_out);
        void getDofForNodes(std::vector<size_t> &Nodes, std::vector<size_t> DofInd, std::vector<PetscInt> &Dofs_out);
        void getDofForNodesSeries(std::vector<size_t> &Nodes, std::vector<size_t> &DofInd, std::vector<size_t> &Dofs_out);
        std::vector<PetscInt> getGhostDofNumbers(size_t Step);

        size_t nDofs;   //number of unique degree of freedom names
        size_t maxSteps;//amount of staggered steps
        std::vector<std::string> DofNames;  //names of degrees of freedom
        std::vector<size_t> DofStep;    //step in which degrees of freedom are resolved

        std::vector<std::vector<size_t>> LocalDofNumRange;  // step, start-stop
        std::vector<size_t> TotalDofRange;
        std::vector<std::vector<std::map<size_t, size_t>>> DofNumbering;  // [step][dofIdx][node]->dofNum;
    private:
        PetscMPIInt    rank, size;

        std::vector<std::vector<PetscInt>> GhostDofs; 
        std::vector<std::vector<size_t>> DofsToAdd;

        std::string FilenamePrefix;
};

#endif