#include "SaveData.h"



/// @brief Initialize saveData file to either load from or write into
/// @param FileName Name of file to create
/// @param OP reading/writing
SaveDataFile::SaveDataFile(std::string FileName, std::string OP){
    PetscInt MPI_size2, MPI_rank2;
    PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&MPI_size2));
    PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&MPI_rank2));
    MPI_size = MPI_size2;
    MPI_rank = MPI_rank2;

    PathPrefix = FileName;
    NamePrefix = "";
    CurrentOP = OP;

    FileAccessProps fapl;
    fapl.add(MPIOFileAccess{PETSC_COMM_WORLD, MPI_INFO_NULL});
    if (CurrentOP == "read"){ //open as read-only
        SaveFile = new File(FileName, File::ReadOnly, fapl);
    } else if (CurrentOP == "write") { //open in read/write mode
        SaveFile = new File(FileName, File::ReadWrite | File::Create | File::Truncate, fapl);
    }
}

SaveDataFile::~SaveDataFile(){
    delete SaveFile;
}

/// @brief Sets a prefix used within the output file
/// @param pfx String indicating the output path to use within the file
void SaveDataFile::SetPrefix(std::string pfx){
    NamePrefix = pfx;
}
/// @brief Gets the current prefix path used within the output file
/// @return Prefix path used within the hdf5 output path
std::string SaveDataFile::GetPrefix(){
    return NamePrefix;
}

/// @brief Saves the provided variable to the save file
/// @param SaveName Name of the to-be-saved variable
/// @param ToSave Pointer to the variable to save
void SaveDataFile::Save(std::string SaveName, double& ToSave){
    DataSet dataset = SaveFile->createDataSet<double>(NamePrefix+"/"+SaveName, DataSpace({MPI_size}));
    dataset.select({MPI_rank},{1}).write(ToSave);
}
/// @brief Loads the provided variable from the save file
/// @param SaveName Name of the to-be-loaded variable
/// @param ToSave Pointer to the variable to overwrite
void SaveDataFile::Load(std::string SaveName, double& ToSave){
    HighFive::DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName);
    dataset.select({MPI_rank},{1}).read(ToSave);
}

/// @brief Saves the provided variable to the save file
/// @param SaveName Name of the to-be-saved variable
/// @param ToSave Pointer to the variable to save
void SaveDataFile::Save(std::string SaveName, bool& ToSave){
    uint tsave = ToSave;
    DataSet dataset = SaveFile->createDataSet<uint>(NamePrefix+"/"+SaveName, DataSpace({MPI_size}));
    dataset.select({MPI_rank},{1}).write(tsave);
}
/// @brief Loads the provided variable from the save file
/// @param SaveName Name of the to-be-loaded variable
/// @param ToSave Pointer to the variable to overwrite
void SaveDataFile::Load(std::string SaveName, bool& ToSave){
    uint tsave;
    DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName);
    dataset.select({MPI_rank},{1}).read(tsave);
    ToSave= tsave;
}

/// @brief Saves the provided variable to the save file
/// @param SaveName Name of the to-be-saved variable
/// @param ToSave Pointer to the variable to save
void SaveDataFile::Save(std::string SaveName, size_t& ToSave){
    DataSet dataset = SaveFile->createDataSet<size_t>(NamePrefix+"/"+SaveName, DataSpace({MPI_size}));
    dataset.select({MPI_rank},{1}).write(ToSave);
}
/// @brief Loads the provided variable from the save file
/// @param SaveName Name of the to-be-loaded variable
/// @param ToSave Pointer to the variable to overwrite
void SaveDataFile::Load(std::string SaveName, size_t& ToSave){
    DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName);
    dataset.select({MPI_rank},{1}).read(ToSave);
}

/// @brief Saves the provided variable to the save file
/// @param SaveName Name of the to-be-saved variable
/// @param ToSave Pointer to the variable to save
void SaveDataFile::Save(std::string SaveName, std::vector<std::vector<size_t>>& ToSave){
    std::vector<size_t> SetRows(MPI_size);
    for (size_t i = 0; i < MPI_size; i++){
        SetRows[i] = 0;
    }
    SetRows[MPI_rank] = ToSave.size();

    MPI_Allreduce(MPI_IN_PLACE, SetRows.data(), MPI_size, my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
    for (size_t i = 0; i < MPI_size; i++){
        std::vector<size_t> SetColumns(SetRows[i]);
        if (i==MPI_rank){
            for (size_t j = 0; j < SetRows[i]; j++){
                SetColumns[j] = ToSave[j].size();
            }
        } else {
            for (size_t j = 0; j < SetRows[i]; j++){
                SetColumns[j] = 0;
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, SetColumns.data(), SetRows[i], my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
        DataSet dataset = SaveFile->createDataSet<size_t>(NamePrefix+"/"+SaveName+"/"+std::to_string(i)+"/nrows", DataSpace(1));
        for (size_t j = 0; j < SetRows[i]; j++){
            DataSet dataset2 = SaveFile->createDataSet<size_t>(NamePrefix+"/"+SaveName+"/"+std::to_string(i)+"/row"+std::to_string(j), DataSpace(SetColumns[j]));
        }
    }
    {
        DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/nrows");
        dataset.write(ToSave.size());
    }
    for (size_t j = 0; j < ToSave.size(); j++){
        DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/row"+std::to_string(j));
        dataset.write(ToSave[j]);
    }
}
/// @brief Loads the provided variable from the save file
/// @param SaveName Name of the to-be-loaded variable
/// @param ToSave Pointer to the variable to overwrite
void SaveDataFile::Load(std::string SaveName, std::vector<std::vector<size_t>>& ToSave){
    size_t nRows;
    SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/nrows").read(nRows);
    ToSave.resize(nRows);
    for (size_t i = 0; i < nRows; i++){
        DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/row"+std::to_string(i));
        dataset.read(ToSave[i]);    
    }
}

/// @brief Saves the provided variable to the save file
/// @param SaveName Name of the to-be-saved variable
/// @param ToSave Pointer to the variable to save
void SaveDataFile::Save(std::string SaveName, std::vector<std::vector<double>>& ToSave){
    std::vector<size_t> SetRows(MPI_size);
    for (size_t i = 0; i < MPI_size; i++){
        SetRows[i] = 0;
    }
    SetRows[MPI_rank] = ToSave.size();

    MPI_Allreduce(MPI_IN_PLACE, SetRows.data(), MPI_size, my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
    for (size_t i = 0; i < MPI_size; i++){
        std::vector<size_t> SetColumns(SetRows[i]);
        if (i==MPI_rank){
            for (size_t j = 0; j < SetRows[i]; j++){
                SetColumns[j] = ToSave[j].size();
            }
        } else {
            for (size_t j = 0; j < SetRows[i]; j++){
                SetColumns[j] = 0;
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, SetColumns.data(), SetRows[i], my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
        DataSet dataset = SaveFile->createDataSet<size_t>(NamePrefix+"/"+SaveName+"/"+std::to_string(i)+"/nrows", DataSpace(1));
        for (size_t j = 0; j < SetRows[i]; j++){
            DataSet dataset2 = SaveFile->createDataSet<double>(NamePrefix+"/"+SaveName+"/"+std::to_string(i)+"/row"+std::to_string(j), DataSpace(SetColumns[j]));
        }
    }
    {
        DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/nrows");
        dataset.write(ToSave.size());
    }
    for (size_t j = 0; j < ToSave.size(); j++){
        DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/row"+std::to_string(j));
        dataset.write(ToSave[j]);
    }
}
/// @brief Loads the provided variable from the save file
/// @param SaveName Name of the to-be-loaded variable
/// @param ToSave Pointer to the variable to overwrite
void SaveDataFile::Load(std::string SaveName, std::vector<std::vector<double>>& ToSave){
    size_t nRows;
    SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/nrows").read(nRows);
    ToSave.resize(nRows);
    for (size_t i = 0; i < nRows; i++){
        DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/row"+std::to_string(i));
        dataset.read(ToSave[i]);    
    }
}

/// @brief Saves the provided variable to the save file
/// @param SaveName Name of the to-be-saved variable
/// @param ToSave Pointer to the variable to save
void SaveDataFile::Save(std::string SaveName, std::vector<std::string>& ToSave){
    std::vector<size_t> SetRows(MPI_size);
    for (size_t i = 0; i < MPI_size; i++){
        SetRows[i] = 0; 
    }
    SetRows[MPI_rank] = ToSave.size();
    MPI_Allreduce(MPI_IN_PLACE, SetRows.data(), MPI_size, my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
    for (size_t i = 0; i < MPI_size; i++){
        DataSet dataset = SaveFile->createDataSet<std::string>(NamePrefix+"/"+SaveName+"/"+std::to_string(i), DataSpace(SetRows[i]));
    }
    DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank));
    dataset.write(ToSave);
}
/// @brief Loads the provided variable from the save file
/// @param SaveName Name of the to-be-loaded variable
/// @param ToSave Pointer to the variable to overwrite
void SaveDataFile::Load(std::string SaveName, std::vector<std::string>& ToSave){
    DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank));
    dataset.read(ToSave);
}

/// @brief Saves the provided variable to the save file
/// @param SaveName Name of the to-be-saved variable
/// @param ToSave Pointer to the variable to save
void SaveDataFile::Save(std::string SaveName, std::vector<size_t>& ToSave){
    std::vector<size_t> SetRows(MPI_size);
    for (size_t i = 0; i < MPI_size; i++){
        SetRows[i] = 0; 
    }
    SetRows[MPI_rank] = ToSave.size();
    MPI_Allreduce(MPI_IN_PLACE, SetRows.data(), MPI_size, my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
    for (size_t i = 0; i < MPI_size; i++){
        DataSet dataset = SaveFile->createDataSet<size_t>(NamePrefix+"/"+SaveName+"/"+std::to_string(i), DataSpace(SetRows[i]));
    }
    DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank));
    dataset.write(ToSave);
}
/// @brief Loads the provided variable from the save file
/// @param SaveName Name of the to-be-loaded variable
/// @param ToSave Pointer to the variable to overwrite
void SaveDataFile::Load(std::string SaveName, std::vector<size_t>& ToSave){
    DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank));
    dataset.read(ToSave);
}

/// @brief Saves the provided variable to the save file
/// @param SaveName Name of the to-be-saved variable
/// @param ToSave Pointer to the variable to save
void SaveDataFile::Save(std::string SaveName, std::vector<std::vector<PetscInt>>& ToSave){
    std::vector<size_t> SetRows(MPI_size);
    for (size_t i = 0; i < MPI_size; i++){
        SetRows[i] = 0;
    }
    SetRows[MPI_rank] = ToSave.size();

    MPI_Allreduce(MPI_IN_PLACE, SetRows.data(), MPI_size, my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
    for (size_t i = 0; i < MPI_size; i++){
        std::vector<size_t> SetColumns(SetRows[i]);
        if (i==MPI_rank){
            for (size_t j = 0; j < SetRows[i]; j++){
                SetColumns[j] = ToSave[j].size();
            }
        } else {
            for (size_t j = 0; j < SetRows[i]; j++){
                SetColumns[j] = 0;
            } 
        }
        MPI_Allreduce(MPI_IN_PLACE, SetColumns.data(), SetRows[i], my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
        DataSet dataset = SaveFile->createDataSet<size_t>(NamePrefix+"/"+SaveName+"/"+std::to_string(i)+"/nrows", DataSpace(1));
        for (size_t j = 0; j < SetRows[i]; j++){
            DataSet dataset = SaveFile->createDataSet<PetscInt>(NamePrefix+"/"+SaveName+"/"+std::to_string(i)+"/row"+std::to_string(j), DataSpace(SetColumns[j]));
        }
    }
    {
        DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/nrows");
        dataset.write(ToSave.size());
    }
    for (size_t j = 0; j < ToSave.size(); j++){
        DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/row"+std::to_string(j));
        dataset.write(ToSave[j]);
    }
}
/// @brief Loads the provided variable from the save file
/// @param SaveName Name of the to-be-loaded variable
/// @param ToSave Pointer to the variable to overwrite
void SaveDataFile::Load(std::string SaveName, std::vector<std::vector<PetscInt>>& ToSave){
    size_t nRows;
    SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/nrows").read(nRows);
    ToSave.resize(nRows);
    for (size_t i = 0; i < nRows; i++){
        DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/row"+std::to_string(i));
        dataset.read(ToSave[i]);    
    }
}

/// @brief Saves the provided variable to the save file
/// @param SaveName Name of the to-be-saved variable
/// @param ToSave Pointer to the variable to save
void SaveDataFile::Save(std::string SaveName, std::vector<std::vector<std::map<size_t, size_t>>>& ToSave){
    std::vector<size_t> SetRows(MPI_size);
    for (size_t i = 0; i < MPI_size; i++){
        SetRows[i] = 0;
    }
    SetRows[MPI_rank] = ToSave.size();

    MPI_Allreduce(MPI_IN_PLACE, SetRows.data(), MPI_size, my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
    for (size_t i = 0; i < MPI_size; i++){
        std::vector<size_t> SetColumns(SetRows[i]);
        if (i==MPI_rank){
            for (size_t j = 0; j < SetRows[i]; j++){
                SetColumns[j] = ToSave[j].size();
            }
        } else {
            for (size_t j = 0; j < SetRows[i]; j++){
                SetColumns[j] = 0;
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, SetColumns.data(), SetRows[i], my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
        DataSet dataset = SaveFile->createDataSet<size_t>(NamePrefix+"/"+SaveName+"/"+std::to_string(i)+"/nPerCol", DataSpace(SetRows[i]));
        for (size_t j = 0; j < SetRows[i]; j++){
            for (size_t k = 0; k < SetColumns[j]; k++){
                std::vector<size_t> MapSize(1);
                if (i == MPI_rank){
                    MapSize[0] = ToSave[j][k].size();
                } else {
                    MapSize[0] = 0;
                }
                MPI_Allreduce(MPI_IN_PLACE, MapSize.data(), 1, my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
                DataSet dataset = SaveFile->createDataSet<size_t>(NamePrefix+"/"+SaveName+"/"+std::to_string(i)+"/map_r"+std::to_string(j)+"_c"+std::to_string(k), DataSpace({MapSize[0],2}));
            }
        }
    }

    std::vector<size_t> SetColumns(ToSave.size());
    for (size_t j = 0; j < ToSave.size(); j++){
        SetColumns[j] = ToSave[j].size();
        for (size_t k = 0; k < ToSave[j].size(); k++){
            std::vector<std::vector<size_t>> MapData(ToSave[j][k].size());
            size_t i=0; 
            for (const auto &it : ToSave[j][k]){
                MapData[i].resize(2);
                MapData[i][0] = it.first;
                MapData[i][1] = it.second;
                i+= 1;
            }   
            if (i>0){
                DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/map_r"+std::to_string(j)+"_c"+std::to_string(k));
                dataset.write(MapData);
            }
        }
    }
    DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/nPerCol");
    dataset.write(SetColumns);
}
/// @brief Loads the provided variable from the save file
/// @param SaveName Name of the to-be-loaded variable
/// @param ToSave Pointer to the variable to overwrite
void SaveDataFile::Load(std::string SaveName, std::vector<std::vector<std::map<size_t, size_t>>>& ToSave){
    std::vector<size_t> SetColumns;
    DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/nPerCol");
    dataset.read(SetColumns);
    ToSave.resize(SetColumns.size());
    for (size_t i = 0; i < SetColumns.size(); i++){
        ToSave[i].resize(SetColumns[i]);
        for (size_t j = 0; j < SetColumns[i]; j++){
            std::vector<std::vector<size_t>> MapData;
            DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/map_r"+std::to_string(i)+"_c"+std::to_string(j));
            dataset.read(MapData);

            for (size_t k = 0; k < MapData.size(); k++){
                ToSave[i][j].insert({MapData[k][0], MapData[k][1]});
            }
        }
    }
}

/// @brief Saves the provided variable to the save file
/// @param SaveName Name of the to-be-saved variable
/// @param ToSave Pointer to the variable to save
void SaveDataFile::Save(std::string SaveName, GhostVector& ToSave){
    std::vector<size_t> GhostSize(MPI_size);
    std::vector<PetscInt> VecSize(MPI_size);
    for (size_t i = 0; i < MPI_size; i++){
        GhostSize[i] = 0;
    }
    GhostSize[MPI_rank] = ToSave.Ghosts.size();
    VecGetLocalSize(ToSave.DataVector, &VecSize[MPI_rank]);

    MPI_Allreduce(MPI_IN_PLACE, GhostSize.data(), MPI_size, my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, VecSize.data(), MPI_size, MPIU_INT, MPI_MAX, PETSC_COMM_WORLD);
    PetscBarrier(NULL);
    for (size_t i = 0; i < MPI_size; i++){
        SaveFile->createDataSet<size_t>(NamePrefix+"/"+SaveName+"/"+std::to_string(i)+"/GhostAlloc", DataSpace(GhostSize[i]));
        SaveFile->createDataSet<size_t>(NamePrefix+"/"+SaveName+"/"+std::to_string(i)+"/LocalSize", DataSpace(1));
        SaveFile->createDataSet<size_t>(NamePrefix+"/"+SaveName+"/"+std::to_string(i)+"/GlobalSize", DataSpace(1));
        SaveFile->createDataSet<double>(NamePrefix+"/"+SaveName+"/"+std::to_string(i)+"/VecData", DataSpace(VecSize[i]));
    }
    
    DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/GhostAlloc");
    dataset.write(ToSave.Ghosts);

    DataSet dataset2 = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/LocalSize");
    PetscInt localSize; VecGetLocalSize(ToSave.DataVector, &localSize);
    dataset2.write(localSize);

    DataSet dataset3 = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/GlobalSize");
    PetscInt globalSize; VecGetSize(ToSave.DataVector, &globalSize);
    dataset3.write(globalSize);

    Vec LocalVec; std::vector<double> VecData(localSize);
    PetscCallThrow(VecGhostGetLocalForm(ToSave.DataVector, &LocalVec));
    PetscScalar* V;
    PetscCallThrow(VecGetArray(LocalVec, &V));
    for (int i = 0; i < localSize; i++){
        VecData[i] = V[i];
    }
    PetscCallThrow(VecRestoreArray(LocalVec, &V));
    PetscCallThrow(VecGhostRestoreLocalForm(ToSave.DataVector,&LocalVec)); 

    DataSet dataset4 = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/VecData");
    dataset4.write(VecData);
}
/// @brief Loads the provided variable from the save file
/// @param SaveName Name of the to-be-loaded variable
/// @param ToSave Pointer to the variable to overwrite
void SaveDataFile::Load(std::string SaveName, GhostVector& ToSave){
    PetscInt LocalSize, GlobalSize;
    std::vector<PetscInt> Ghosts;
    DataSet dataset = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/GhostAlloc");
    dataset.read(Ghosts);

    DataSet dataset2 = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/LocalSize");
    dataset2.read(LocalSize);

    DataSet dataset3 = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/GlobalSize");
    dataset3.read(GlobalSize);

    ToSave.SetSize(LocalSize, GlobalSize, Ghosts);

    std::vector<double> VecData;
    DataSet dataset4 = SaveFile->getDataSet(NamePrefix+"/"+SaveName+"/"+std::to_string(MPI_rank)+"/VecData");
    dataset4.read(VecData);

    Vec LocalVec;
    PetscCallThrow(VecGhostGetLocalForm(ToSave.DataVector, &LocalVec));
    PetscScalar* V;
    PetscCallThrow(VecGetArray(LocalVec, &V));
    for (int i = 0; i < LocalSize; i++){
        V[i] = VecData[i];
    }
    PetscCallThrow(VecRestoreArray(LocalVec, &V));
    PetscCallThrow(VecGhostRestoreLocalForm(ToSave.DataVector,&LocalVec)); 

    ToSave.AssemblyStart();
    ToSave.AssemblyEnd();
    ToSave.SyncStart(INSERT_VALUES);
    ToSave.SyncEnd(INSERT_VALUES);
}

