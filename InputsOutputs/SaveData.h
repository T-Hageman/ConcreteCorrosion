#ifndef SAVEDATA_H
#define SAVEDATA_H

#include <petsc.h>
#include <fstream>
#include <iostream>
#include <vector>

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

using namespace HighFive;

#include "inputData.h"
#include "../utility/utility.h"

/// @brief Class to handle parrallel writing of restart data to a single file
class SaveDataFile {
    public:
        SaveDataFile(std::string FileName, std::string OP);
        ~SaveDataFile();

        void SetPrefix(std::string pfx);
        std::string GetPrefix();

        void Save(std::string SaveName, double& ToSave);
        void Load(std::string SaveName, double& ToSave);

        void Save(std::string SaveName, bool& ToSave);
        void Load(std::string SaveName, bool& ToSave);

        void Save(std::string SaveName, size_t& ToSave);
        void Load(std::string SaveName, size_t& ToSave);

        void Save(std::string SaveName, std::vector<std::vector<size_t>>& ToSave);
        void Load(std::string SaveName, std::vector<std::vector<size_t>>& ToSave);

        void Save(std::string SaveName, std::vector<std::vector<double>>& ToSave);
        void Load(std::string SaveName, std::vector<std::vector<double>>& ToSave);

        void Save(std::string SaveName, std::vector<size_t>& ToSave);
        void Load(std::string SaveName, std::vector<size_t>& ToSave);

        void Save(std::string SaveName, std::vector<std::string>& ToSave);
        void Load(std::string SaveName, std::vector<std::string>& ToSave);

        void Save(std::string SaveName, std::vector<std::vector<PetscInt>>& ToSave);
        void Load(std::string SaveName, std::vector<std::vector<PetscInt>>& ToSave);

        void Save(std::string SaveName, std::vector<std::vector<std::map<size_t, size_t>>>& ToSave);
        void Load(std::string SaveName, std::vector<std::vector<std::map<size_t, size_t>>>& ToSave);

        void Save(std::string SaveName, GhostVector& ToSave);
        void Load(std::string SaveName, GhostVector& ToSave);
    private:
        std::string NamePrefix; //Internal file prefix which gets prepended to the savename
        std::string PathPrefix; //path and savefile name prefix

        size_t MPI_rank, MPI_size;

        std::string CurrentOP;  //read/write
        File* SaveFile;         //pointer to hdf savefile
};

#endif