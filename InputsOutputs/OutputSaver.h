#ifndef OUTPUTSAVER_H
#define OUTPUTSAVER_H

#include <petsc.h>
#include <fstream>
#include <iostream>

#include "../Physics/physics.h"
#include "inputData.h"

/// @brief Object to save outputs to a .hdf5 file in parrallel
class OutputDataSaver {
    public:
        OutputDataSaver(Physics& physics, inputData& inputs);
        ~OutputDataSaver();

        void saveResults(size_t timeStep, double time);
        void saveMesh(size_t timeStep);

    private:
        bool OutputMeshOnlyOnce; // Flag indicating whether to save the mesh once to a separate file (true), or include within each output file (false)
        std::vector<std::vector<size_t>> MyRange;   //Owned range within the current core, for which it should save outputs
        std::vector<size_t> TotalSize;  //Total size of each output group (=number of elements within group)

        std::string FilenamePrefix; //prefix or filePath added to saved output files
        std::vector<std::string> GroupNamesToSave; //array of element group names for which outputs are provided
        std::vector<size_t> GroupIndicesToSave; //array of element group indices for which outputs are provided
        std::vector<std::vector<std::string>> NodalQuantityToSave, IPQuantitiesToSave; // arrays of names of quantities to save [Group][quantityIndex]->quantityName

        Physics* physics;   //pointer to physics object
        Mesh* mesh;         //pointer to mesh
        PetscMPIInt    MPIsize, MPIrank;
};

#endif