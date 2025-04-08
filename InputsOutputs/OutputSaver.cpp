#include "OutputSaver.h"

#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5File.hpp>

using namespace HighFive;

/// @brief Initialization, sets up sizes and counters for parrallel saving
/// @param physics_in Pointer to physics object
/// @param inputs Input data object
OutputDataSaver::OutputDataSaver(Physics& physics_in, inputData& inputs){
    PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&MPIsize));
    PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&MPIrank));

    physics = &physics_in;
    mesh = physics->mesh;

    inputs.GetRequired(FilenamePrefix,{"Outputs","SaveFolder"});
    inputs.GetRequired(OutputMeshOnlyOnce,{"Outputs","OutputMeshOnlyOnce"});

    // Get the groups and quantaties to save for these groups
    inputs.GetRequired(GroupNamesToSave, {"Outputs","ElementGroups"});

    GroupIndicesToSave.resize(GroupNamesToSave.size());
    NodalQuantityToSave.resize(GroupNamesToSave.size());
    IPQuantitiesToSave.resize(GroupNamesToSave.size());

    // Get range contained within each core
    MyRange.resize(GroupNamesToSave.size());
    TotalSize.resize(GroupNamesToSave.size());
    for (size_t i = 0; i < GroupIndicesToSave.size(); i++){
        GroupIndicesToSave[i] = physics->mesh->GetElementGroupIdx(GroupNamesToSave[i]);

        inputs.GetRequired(NodalQuantityToSave[i], {"Outputs",GroupNamesToSave[i],"Nodes"});
        inputs.GetRequired(IPQuantitiesToSave[i], {"Outputs",GroupNamesToSave[i],"IP"});

        MyRange[i].resize(2);
        std::vector<size_t> Offset(MPIsize); for (PetscMPIInt j = 0; j < MPIsize; j++) Offset[j] = 0;
        Offset[MPIrank] = mesh->ElementGroups[GroupIndicesToSave[i]].NElems;
        MPI_Allreduce(MPI_IN_PLACE, Offset.data(), MPIsize, my_MPI_SIZE_T, MPI_SUM, PETSC_COMM_WORLD);
        TotalSize[i] = 0;
        MyRange[i][0] = 0;
        for (PetscMPIInt j = 0; j < MPIrank; j++){
            MyRange[i][0] += Offset[j];
        }
        MyRange[i][1] = MyRange[i][0] + Offset[MPIrank];
        //std::cout << MyRange[i][0] << "   " << MyRange[i][1] << "\n";
        for (PetscMPIInt j = 0; j < MPIsize; j++){
            TotalSize[i] += Offset[j];
        }
    }

    saveMesh(0);
}

OutputDataSaver::~OutputDataSaver(){

}

void OutputDataSaver::saveResults(size_t timeStep, double time){
    Logs.PrintSingle("Outputting results to File:\n",1);
    if (OutputMeshOnlyOnce == false){
        saveMesh(timeStep);
    }

    std::string FileName = FilenamePrefix + "/results_" + std::to_string(timeStep) + ".hdf5";

    FileAccessProps fapl;
    fapl.add(MPIOFileAccess{PETSC_COMM_WORLD, MPI_INFO_NULL});

    File file(FileName, File::ReadWrite | File::Create | File::Truncate, fapl);

    //Indicidual data to save
    DataSet dataset = file.createDataSet<double>("time", DataSpace(1));
    if (MPIrank == 0){
        dataset.write(time);
    }

    //Distributed data to save
    for (size_t eg=0; eg<GroupIndicesToSave.size(); eg++){
        size_t egInd = GroupIndicesToSave[eg];

        size_t nNodes = mesh->ElementGroups[egInd].BaseElem->NExport.size();
        size_t nip = mesh->ElementGroups[egInd].BaseElem->ipcount;

        std::vector<std::vector<double>> NodeData, IPData;
        NodeData.resize(mesh->ElementGroups[egInd].NElems); 
        for (size_t i = 0; i < mesh->ElementGroups[egInd].NElems; i++) NodeData[i].resize(nNodes);
        IPData.resize(mesh->ElementGroups[egInd].NElems); 
        for (size_t i = 0; i < mesh->ElementGroups[egInd].NElems; i++) IPData[i].resize(nip);

        for (size_t i = 0; i < NodalQuantityToSave[eg].size(); i++){
            physics->GetNodalDataToSave(NodalQuantityToSave[eg][i], egInd, NodeData);
            DataSet dataset = file.createDataSet<double>(GroupNamesToSave[eg]+"/"+NodalQuantityToSave[eg][i], DataSpace({TotalSize[eg], nNodes}));
            if (MyRange[eg][1]-MyRange[eg][0]>0){
                dataset.select({MyRange[eg][0],0},{MyRange[eg][1]-MyRange[eg][0],nNodes} ).write(NodeData);
            }
        }
        for (size_t i = 0; i < IPQuantitiesToSave[eg].size(); i++){
            physics->GetIPDataToSave(IPQuantitiesToSave[eg][i], egInd, IPData);
            DataSet dataset = file.createDataSet<double>(GroupNamesToSave[eg]+"/"+IPQuantitiesToSave[eg][i], DataSpace({TotalSize[eg], nip}));
            if (MyRange[eg][1]-MyRange[eg][0]>0){
                dataset.select({MyRange[eg][0],0},{MyRange[eg][1]-MyRange[eg][0],nip} ).write(IPData);
            }
        }
    }
    Logs.PrintSingle("Finished outputting results to File\n",1);
 
    if (MPIrank == 0){
        Logs.PrintSingle("Outputting timedata to File:\n",1);

        std::string FileName = FilenamePrefix + "/TimeData.hdf5";
        File fileTimeSeries(FileName, File::ReadWrite | File::Create | File::Truncate );
        size_t nTimeData = physics->TimeDataTypes.size();
        size_t lTimeData = physics->TimeData[0].size();

        DataSet datasetTimeTypes  = fileTimeSeries.createDataSet<std::string>("TimeDataTypes", DataSpace(nTimeData));
        DataSet datasetTimeSeries = fileTimeSeries.createDataSet<double>("TimeData", DataSpace(nTimeData, lTimeData));

        datasetTimeTypes.write(physics->TimeDataTypes);
        datasetTimeSeries.write(physics->TimeData);

        Logs.PrintSingle("Finished outputting timedata to File\n",1);
    }
    
}

/// @brief Saves the mesh locations in an easy-to-use file format for post-processing
/// @param timeStep Timestep number to append to output file
void OutputDataSaver::saveMesh(size_t timeStep){
    Logs.PrintSingle("Outputting mesh to File:\n",1);

    //open file
    std::string FileName = FilenamePrefix + "/mesh_" + std::to_string(timeStep)+".hdf5";

    FileAccessProps fapl;
    fapl.add(MPIOFileAccess{PETSC_COMM_WORLD, MPI_INFO_NULL});

    File file(FileName, File::ReadWrite | File::Create | File::Truncate, fapl);

    //loop over output groups
    for (size_t eg=0; eg<GroupIndicesToSave.size(); eg++){  
        size_t egInd = GroupIndicesToSave[eg];
        size_t nNodes, nip;
        nip = mesh->ElementGroups[egInd].BaseElem->ipcount;
        std::vector<std::vector<double>> XCoords, YCoords, ZCoords, XIPCoords, YIPCoords, ZIPCoords;

        nNodes = mesh->GetExportMeshCoords(egInd, XCoords, YCoords, ZCoords, XIPCoords, YIPCoords, ZIPCoords);

        //save to output file
        DataSet dataset_X = file.createDataSet<double>(GroupNamesToSave[eg]+"/X", DataSpace({TotalSize[eg], nNodes}));
        DataSet dataset_Y = file.createDataSet<double>(GroupNamesToSave[eg]+"/Y", DataSpace({TotalSize[eg], nNodes}));
        DataSet dataset_XIP = file.createDataSet<double>(GroupNamesToSave[eg]+"/Xip", DataSpace({TotalSize[eg], nip}));
        DataSet dataset_YIP = file.createDataSet<double>(GroupNamesToSave[eg]+"/Yip", DataSpace({TotalSize[eg], nip}));

        if (MyRange[eg][1]-MyRange[eg][0]>0){
            dataset_X.select({MyRange[eg][0],0},{MyRange[eg][1]-MyRange[eg][0],nNodes} ).write(XCoords);
            dataset_Y.select({MyRange[eg][0],0},{MyRange[eg][1]-MyRange[eg][0],nNodes} ).write(YCoords);
            dataset_XIP.select({MyRange[eg][0],0},{MyRange[eg][1]-MyRange[eg][0],nip} ).write(XIPCoords);
            dataset_YIP.select({MyRange[eg][0],0},{MyRange[eg][1]-MyRange[eg][0],nip} ).write(YIPCoords);
        }

        if (mesh->dim==3){
            DataSet dataset_Z = file.createDataSet<double>(GroupNamesToSave[eg]+"/Z", DataSpace({TotalSize[eg], nNodes}));
            DataSet dataset_ZIP = file.createDataSet<double>(GroupNamesToSave[eg]+"/Zip", DataSpace({TotalSize[eg], nip}));  
            if (MyRange[eg][1]-MyRange[eg][0]>0){
                dataset_Z.select({MyRange[eg][0],0},{MyRange[eg][1]-MyRange[eg][0],nNodes} ).write(ZCoords);
                dataset_ZIP.select({MyRange[eg][0],0},{MyRange[eg][1]-MyRange[eg][0],nip} ).write(ZIPCoords);
            }
        }
    }
    Logs.PrintSingle("Finished outputting mesh to File\n",1);
}




