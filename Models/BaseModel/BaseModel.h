#ifndef BASEMODEL_H
#define BASEMODEL_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../../InputsOutputs/inputData.h"
#include "../../mesh/mesh.h"
#include "../../Physics/DofSpace.h"
#include "../../Physics/Constrainer.h"
#include "../../InputsOutputs/SaveData.h"
#include "../../utility/utility.h"

namespace CommitTypes{
    constexpr int TIMEDEP_COMMIT_TYPE = 1;
    constexpr int PATHDEP_COMMIT_TYPE = 2;
}

class Physics;

/// @brief This class implements the interfaces used by the physics models
/** 
Within the initialization of simulations, all physics models are added based on
the name used. For this inheritance to work for any named model, they need to
inherit this BaseModel class, and follow/overwrite the required functions.
Please note that within derived models, not all functions need to be declared,
leaving blank functions just skips them (e.g. when a model does not add any time
data to save). 

*/ 
class BaseModel {
    public:
        std::string MyName;     // specific name that of this model (as it is called in the input file)
        std::string ModelName;  // general name of this model
        uint dim;

        /// @brief Creates the model and links it to the overall physics
        /// objects (physics, mesh, dofspace). 
        /// @param My_Physics Physics object which will be linked to this model
        /// @param MyName Name used within input files to designate this model
        BaseModel(Physics& My_Physics, std::string MyName);

        virtual ~BaseModel();

        /// @brief Reloads this model from a previous restart file. Please note
        /// that most parameters do not to be saved, and can be re-initialized
        /// from the input file. Only history-type variables need to be saved/loaded
        /// @param inputs input file where properties are taken from
        /// @param data Save file from which model is being restored
        virtual void load(inputData& inputs, SaveDataFile& data);

        /// @brief Saves variables of this model to a restart file, for later
        /// restarting/error recovery
        /// @param data Data object where to save the restart data to.
        virtual void save(SaveDataFile& data);

        /// @brief Initializes a model for the first time. This function is only
        /// ever called once per model, with re-starts calling load instead.
        /// This function is run after the object is initialized via the constructor.
        /// @param inputs Input json file structure.
        virtual void init(inputData& inputs);

        /// @brief Request that this model adds the relevant physics to the
        /// tangent stiffness matrix K and force vector f for the current step.
        /// This function always gets called, independent of whtehr the model
        /// adds anything during the current solving step. As such, the model
        /// should check step against the steps in which their relevant physics
        /// are solved
        /// @param K Tangent stiffness matrix of step
        /// @param f Force vector of step
        /// @param cons constraints that are being applied in the current step
        /// @param step step that is currently being solved
        virtual void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);

        /// @brief Function that gets called once the current step is converged.
        /// Use this to move internal history variables forward in time. 
        /// @param CommitType Reason why the history variables should be
        /// committed, either TIMEDEP_COMMIT_TYPE for the end of a time
        /// increment, or PATHDEP_COMMIT_TYPE if the step is still ongoing, but
        /// other physics perform actions which could cause loading and
        /// subsequent unloading within a single step
        virtual void Commit(int CommitType);

        /// @brief Function To reset the progress made towards the current time/load increment
        virtual void ResetStep();

        /// @brief Saves output data to a savefile for post-processing. Degrees
        /// of freedom are auto-saved, but this function allows derived
        /// parameters to also be saved (e.g. stresses/fluxes).
        /// @param SaveLoc either "Nodes" or "ip" to indicate nodal data or
        /// interpolation point level data should be saved
        /// @param DataName Name of what is being asked to save
        /// @param ElemGroup Name of the element group that should be saved for,
        /// check that this corresponds to the group on which the model operates
        /// @param Data object in which to save the data to, in the format data[element][point]
        /// @return Has any data been added to data
        virtual bool SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data);

        /// @brief Check if this model keeps any time-series which should be outputted
        /// @param DataNames Resize this to add any time data that will be saved
        /// @return number of time data doubles being added
        virtual size_t hasTimeData(std::vector<std::string>& DataNames);

        /// @brief Saves any time-series data of the last time step
        /// @param DataValues Values to be saved
        virtual void GetTimeData(std::vector<double>& DataValues);

    protected:
        Physics* physics;   // pointer to the physics objects that containes state-vectors etc.
        Mesh* mesh;         // pointer to the mesh object
        DofSpace* dofs;     // pointer to the object containing degree of freedom numberings

        PetscMPIInt MPI_rank, MPI_size;

    private:

};

void Register_BaseModel();
BaseModel* New_BaseModel(Physics& My_Physics, std::string MyNameIn);

#endif

