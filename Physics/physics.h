#ifndef PHYSICS_H
#define PHYSICS_H

#include <vector>
#include <Eigen/Dense>
#include <petscvec.h>

#include "../mesh/mesh.h"
#include "../utility/utility.h"
#include "../InputsOutputs/inputData.h"
#include "../Models/ModelRegister.h"
#include "DofSpace.h"
#include "Constrainer.h"
#include "TimeSchemes.h"

#ifdef LOCALENVIRONMENT
    #include "Visualisation/Visualisation.h"
#endif

/// @brief Class which handles the physics descriptions of the problem (e.g. assembling force and stiffness matrix, etc.)
class Physics {
    public:
        Physics(inputData& inputs, Mesh& inMesh);
        Physics(inputData& inputs, Mesh& inMesh, SaveDataFile& data);
        void Save(SaveDataFile& data);
        ~Physics();

        Mesh* mesh; //pointer to mesh
        DofSpace* dofspace; //pointer to degrees of freedom

        size_t nModels; //number of physics models
        std::vector<std::string> ModelNames, ModelTypeNames;    //unique and characteristic names of physics models
        std::vector<BaseModel*> Models; //array of initialized physics models
        std::vector<GhostVector> StateVectors, StateVectorsOld, dStateVectors, dStateVectorsOld, ddStateVectors, ddStateVectorsOld; //system state vectors

        void InitKf();
        void Assemble(size_t step);
        void Constraint(size_t step);
        void AssembleFOnly(size_t step);
        void ConstraintFOnly(size_t step);

        void UpdateState(size_t step, Vec& dx);
        void ResetState();
        void UpdateVelAcc();
        void Commit(int CommitType); // TIMEDEP_COMMIT_TYPE = 1   or   PATHDEP_COMMIT_TYPE = 2;
        void BoundaryConsToOld();

        std::vector<Vec> f, fCon;   //Force vectors for each staggered step
        std::vector<Mat> K, KCon;   //stiffness matrix for each staggered step
        std::vector<Vec> CForceVec, CForceVec2; //Constrained forces
        Constrainer* cons;  //constraints handling
        TimeScheme* timeScheme; //Time discretisation schemes
        double time;    //current time
        double MaxTime;

        void Init_Timeseries();
        void Append_Timeseries();
        double GetTimeData(std::string DataName);
        std::vector<std::string> TimeDataTypes;
        std::vector<std::vector<double>> TimeData;

        void GetNodalDataToSave(std::string NodalQuantityToSave, size_t EGroupIndex, std::vector<std::vector<double>>& NodeData);
        void GetIPDataToSave(std::string IPQuantityToSave, size_t EGroupIndex, std::vector<std::vector<double>>& IPData);

        #ifdef LOCALENVIRONMENT
            Visualisation* Vis_Plots;
        #endif
    private:
        PetscMPIInt    rank, size;

        void InitStateVecs();
        size_t nAlloc = 1500;    //number of non-zero values per row of the stiffness matrix (initial guess)

        std::vector<size_t> dofVersionUsed; //flag to check whether to re-update constraints

        bool exportMatrices = false; //debugging, export system matrices to file
};

#endif