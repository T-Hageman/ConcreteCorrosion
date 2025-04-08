#include "DofSpace.h"
#include "../utility/utility.h"

#include <petsc.h>
#include <iostream>
#include <string>
#include <highfive/H5File.hpp>

DofSpace::DofSpace(inputData& inputs) {
    PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));

    inputs.GetRequired(DofNames, {"Dofs","DofNames"});
    inputs.GetRequired(DofStep, {"Dofs","DofStep"});
    nDofs = DofNames.size();
    if (nDofs!=DofStep.size()){
        throw std::invalid_argument("Please specify dof steps for all dofs");
    }

    maxSteps = *std::max_element(DofStep.begin(), DofStep.end())+1;
    DofNumbering.resize(maxSteps);
    DofsToAdd.resize(nDofs);
    GhostDofs.resize(maxSteps);
    LocalDofNumRange.resize(maxSteps);
    TotalDofRange.resize(maxSteps);

    inputs.GetRequired(FilenamePrefix,{"Outputs","SaveFolder"});
}

DofSpace::~DofSpace() {

}

/// @brief Saves dofspace to restartable data file
/// @param data 
void DofSpace::save(SaveDataFile& data){
    data.Save("LocalDofNumRange", LocalDofNumRange);
    data.Save("TotalDofRange", TotalDofRange);
    data.Save("DofNumbering", DofNumbering);
    data.Save("GhostDofs", GhostDofs);
}

/// @brief Load dofspace from restartable data file
/// @param data 
void DofSpace::load(SaveDataFile& data){
    data.Load("LocalDofNumRange", LocalDofNumRange);
    data.Load("TotalDofRange", TotalDofRange);
    data.Load("DofNumbering", DofNumbering);
    data.Load("GhostDofs", GhostDofs);
}

/// @brief For a provided dof name, return the index and step in which they are resolved
/// @param names input: Strings to identify the type of degree of freedom
/// @param DofTypes_out output: Index of the degree of freedom type
/// @param DofSteps_out output: Step in which degree of freedom is resolved
void DofSpace::getDofTypesSteps(std::vector<std::string> names, std::vector<size_t> &DofTypes_out, std::vector<size_t> &DofSteps_out){
    DofTypes_out.resize(names.size());
    DofSteps_out.resize(names.size());

    for (size_t i = 0; i < names.size(); i++){
        getDofTypesSteps(names[i], DofTypes_out[i], DofSteps_out[i]);
    }
    return;
}

/// @brief For a provided dof name, return the index and step in which they are resolved
/// @param name input: Strings to identify the type of degree of freedom
/// @param DofType_out output: Index of the degree of freedom type
/// @param DofStep_out output: Step in which degree of freedom is resolved
void DofSpace::getDofTypesSteps(std::string name, size_t &DofType_out, size_t &DofStep_out){
    bool HasDOF = hasDofType(name, DofType_out, DofStep_out);
    if (HasDOF == false){
        throw std::invalid_argument("Dof of type \"" + name + "\" is not defined");
    }  
}

/// @brief Check if degree of freedom exists, if it does, provide info about it
/// @param name input: Strings to identify the type of degree of freedom
/// @param DofType_out output: Index of the degree of freedom type (only defined if output = true)
/// @param DofStep_out output: Step in which degree of freedom is resolved (only defined if output = true)
/// @return does this dof exist, yes/no
bool DofSpace::hasDofType(std::string name, size_t &DofType_out, size_t &DofStep_out){
    auto it = std::find(DofNames.begin(), DofNames.end(), name);
    if (it != DofNames.end()) {
        DofType_out = it - DofNames.begin();
        DofStep_out = DofStep[DofType_out];
        return true;
    } else {
        return false;
    }    
}

/// @brief Check if degree of freedom exists
/// @param name input: Strings to identify the type of degree of freedom
/// @return does this dof exist, yes/no
bool DofSpace::hasDofType(std::string name){
    size_t dummy1, dummy2;
    return hasDofType(name, dummy1, dummy2);
}

/// @brief Selects degrees of freedom to be added to nodes
/// @param Nodes Nodes to which dofs should be added
/// @param DofInd Index of the dof to add
void DofSpace::AddDofs(std::vector<size_t> &Nodes, std::vector<size_t> &DofInd){
    for (size_t idof = 0; idof < DofInd.size(); idof++){
        AddDofs(Nodes, DofInd[idof]);
    }
}

/// @brief Selects degrees of freedom to be added to nodes
/// @param Nodes Nodes to which dofs should be added
/// @param DofInd Index of the dof to add
void DofSpace::AddDofs(std::vector<size_t> &Nodes, size_t DofInd){
    DofsToAdd[DofInd].insert(DofsToAdd[DofInd].end(), Nodes.begin(), Nodes.end());
    std::sort(DofsToAdd[DofInd].begin(), DofsToAdd[DofInd].end());
    auto last = std::unique(DofsToAdd[DofInd].begin(), DofsToAdd[DofInd].end());
    DofsToAdd[DofInd].erase(last, DofsToAdd[DofInd].end());
}

/// @brief Selects degrees of freedom to be added to nodes
/// @param Node Nodes to which dofs should be added
/// @param DofInd Index of the dof to add
void DofSpace::AddDofs(size_t Node, size_t DofInd){
    DofsToAdd[DofInd].insert(DofsToAdd[DofInd].end(), Node);
    std::sort(DofsToAdd[DofInd].begin(), DofsToAdd[DofInd].end());
    auto last = std::unique(DofsToAdd[DofInd].begin(), DofsToAdd[DofInd].end());
    DofsToAdd[DofInd].erase(last, DofsToAdd[DofInd].end());
}

/// @brief Flushes degrees of freedom to a permanent structure, and syncs across nodes
/// @param mesh pointer to mesh (to get owner ranges)
void DofSpace::SyncDofs(Mesh* mesh){
    // Get Node Owner Ranges
    std::vector<size_t> NodeRanges(size+1);
    const PetscInt *NodeRangesPETSC;
    PetscCallThrow(VecGetOwnershipRanges(mesh->Xcoords.DataVector, &NodeRangesPETSC));
    for (int i = 0; i < size+1; i++){
        NodeRanges[i] = NodeRangesPETSC[i];
    }

    //send requested dofs to other cores to verify they are all added (doing this first to not have to fix dof numbering later on)
    {
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> RequestNum(size,size); RequestNum.setZero();
        Eigen::Matrix<std::vector<size_t>, Eigen::Dynamic, Eigen::Dynamic> RequestNodes(size,size), RequestDofs(size,size);
        for (size_t dofInd = 0; dofInd < nDofs; dofInd++){
            for (size_t nd = 0; nd < DofsToAdd[dofInd].size(); nd++){
                int c_to_add_to = 0;
                for (int c = 0; c < size; c++){
                    if (DofsToAdd[dofInd][nd]>=NodeRanges[c]){
                        c_to_add_to = c;
                    }
                }
                if (c_to_add_to!=rank){
                    RequestNum(rank,c_to_add_to) += 1;
                    RequestNodes(rank,c_to_add_to).push_back(DofsToAdd[dofInd][nd]);
                    RequestDofs(rank,c_to_add_to).push_back(dofInd);
                }
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, RequestNum.data(), size*size, my_MPI_SIZE_T, MPI_MAX, PETSC_COMM_WORLD);

        for (size_t c = 0; c < size; c++){
            MPI_Status mySendRequests;
            if (c==rank){ //receiving and adding
                for (size_t c2 = 0; c2 < size; c2++){
                    if (c2!=rank){
                        size_t nToCheck = RequestNum(c2, rank);
                        std::vector<size_t> RecNodes(nToCheck), RecDofs(nToCheck);
                        MPI_Recv(RecNodes.data(), nToCheck, my_MPI_SIZE_T, c2, c2*10+1, PETSC_COMM_WORLD, &mySendRequests);
                        MPI_Recv(RecDofs.data(), nToCheck, my_MPI_SIZE_T, c2, c2*10+2, PETSC_COMM_WORLD, &mySendRequests);

                        for (size_t i = 0; i < RecNodes.size(); i++){
                            AddDofs(RecNodes[i], RecDofs[i]);
                        }
                    }
                }
            } else { //sending
                size_t n_to_send = RequestNum(rank, c);
                MPI_Send(RequestNodes(rank, c).data(), n_to_send, my_MPI_SIZE_T, c, rank*10+1, PETSC_COMM_WORLD);
                MPI_Send(RequestDofs(rank, c).data(), RequestNum(rank, c), my_MPI_SIZE_T, c, rank*10+2, PETSC_COMM_WORLD);
            }
        }
    }

    //set up dof space
    for (size_t CurdofStep = 0; CurdofStep < maxSteps; CurdofStep++){
        //sort Dofs to Add into "buckets"
        std::vector<std::vector<std::vector<size_t>>> DofsToAddPerCore; DofsToAddPerCore.resize(size);   //(Core;dof -> Nodenumbers)
        for (int i = 0; i < size; i++){
            DofsToAddPerCore[i].resize(nDofs);
        }
        
        Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> RequestNum(size,size); RequestNum.setZero();
        Eigen::Matrix<std::vector<size_t>, Eigen::Dynamic, Eigen::Dynamic> RequestNodes(size,size), RequestDofs(size,size);
        for (int i = 0; i < size; i++){
            for (int j = 0; j < size; j++){
                RequestNodes(i,j).resize(0);
            } 
        }
        
        for (size_t dofInd = 0; dofInd < nDofs; dofInd++){
            if (DofStep[dofInd] == CurdofStep){
                for (size_t nd = 0; nd < DofsToAdd[dofInd].size(); nd++){
                    int c_to_add_to = 0;
                    for (int c = 0; c < size; c++){
                        if (DofsToAdd[dofInd][nd]>=NodeRanges[c]){
                            c_to_add_to = c;
                        }
                    }
                    DofsToAddPerCore[c_to_add_to][dofInd].insert(DofsToAddPerCore[c_to_add_to][dofInd].end(), DofsToAdd[dofInd][nd]);
                    if (c_to_add_to!=rank){
                        RequestNum(rank,c_to_add_to) += 1;
                        RequestNodes(rank,c_to_add_to).push_back(DofsToAdd[dofInd][nd]);
                        RequestDofs(rank,c_to_add_to).push_back(dofInd);
                    }
                }
            }
        }

        std::vector<size_t> TotalDofs(size); for (int i = 0; i < size; i++) TotalDofs[i] = 0;
        for (int i = 0; i < size; i++){
            for (size_t j = 0; j < nDofs; j++){
                TotalDofs[i] += DofsToAddPerCore[i][j].size();
            }
        }

        // sync dofs per core and request table
        std::vector<size_t> OwnedDofs(size); for (int i = 0; i < size; i++) OwnedDofs[i] = 0;
        OwnedDofs[rank] = TotalDofs[rank];
        MPI_Allreduce(MPI_IN_PLACE, OwnedDofs.data(), size, my_MPI_SIZE_T, MPI_SUM, PETSC_COMM_WORLD);
        MPI_Allreduce(MPI_IN_PLACE, RequestNum.data(), size*size, my_MPI_SIZE_T, MPI_SUM, PETSC_COMM_WORLD);
        PetscCallThrow(PetscBarrier(NULL)); 

        //set up local dof numberings
        DofNumbering[CurdofStep].resize(nDofs);

        LocalDofNumRange[CurdofStep].resize(2);  LocalDofNumRange[CurdofStep][0] = 0; LocalDofNumRange[CurdofStep][1] = 0;
        TotalDofRange[CurdofStep] = 0;
        for (int i = 0; i < rank; i++){
            LocalDofNumRange[CurdofStep][0] += OwnedDofs[i];
        }
        for (int i = 0; i < rank+1; i++){
            LocalDofNumRange[CurdofStep][1] += OwnedDofs[i];
            
        }
        for (int i = 0; i < size; i++){
            TotalDofRange[CurdofStep] += OwnedDofs[i];
        }
        
        size_t DofCounter = LocalDofNumRange[CurdofStep][0];
        for (size_t nDof = 0; nDof < nDofs; nDof++){
            for (size_t dof = 0; dof < DofsToAddPerCore[rank][nDof].size(); dof++){
                DofNumbering[CurdofStep][nDof].insert({DofsToAddPerCore[rank][nDof][dof],DofCounter});
                DofCounter += 1;
            }
        }
        
        if (size>1){// Get ghosted dof numbering
            //send out requests
            std::vector<MPI_Request> mySendRequests; mySendRequests.reserve(2*size*size+10); 
            std::vector<size_t> mySendIDs; mySendIDs.reserve(2*size*size+10); 
            std::vector<size_t> myReceiveIDs; myReceiveIDs.reserve(2*size*size+10); 
            for (int RequestTo = 0; RequestTo < size; RequestTo++){
                if (RequestNum(rank, RequestTo)>0){
                    mySendRequests.push_back(0); mySendIDs.push_back(rank*size+RequestTo);
                    MPI_Isend(RequestNodes(rank,RequestTo).data(), RequestNodes(rank,RequestTo).size(), my_MPI_SIZE_T, RequestTo, mySendIDs.back(), PETSC_COMM_WORLD, &mySendRequests.back());
                    mySendRequests.push_back(0); mySendIDs.push_back(size*size+rank*size+RequestTo);
                    MPI_Isend(RequestDofs(rank,RequestTo).data(), RequestDofs(rank,RequestTo).size(), my_MPI_SIZE_T, RequestTo, mySendIDs.back(), PETSC_COMM_WORLD, &mySendRequests.back());
                }
            }
            //receive requests
            for (int RequestFrom = 0; RequestFrom < size; RequestFrom++){
                if (RequestNum(RequestFrom, rank)>0){
                    RequestNodes(RequestFrom,rank).resize(RequestNum(RequestFrom, rank));
                    RequestDofs(RequestFrom,rank).resize(RequestNum(RequestFrom, rank));

                    mySendRequests.push_back(0);  myReceiveIDs.push_back(RequestFrom*size+rank);
                    MPI_Irecv(RequestNodes(RequestFrom,rank).data(), RequestNum(RequestFrom, rank), my_MPI_SIZE_T, RequestFrom, myReceiveIDs.back(), PETSC_COMM_WORLD, &mySendRequests.back());
                    mySendRequests.push_back(0);  myReceiveIDs.push_back(size*size+RequestFrom*size+rank);
                    MPI_Irecv(RequestDofs(RequestFrom,rank).data(), RequestNum(RequestFrom, rank), my_MPI_SIZE_T, RequestFrom, myReceiveIDs.back(), PETSC_COMM_WORLD, &mySendRequests.back());
                }
            }
            std::vector<MPI_Status> myStatus(mySendRequests.size());
            MPI_Waitall(mySendRequests.size(), mySendRequests.data(), myStatus.data());
            PetscCallThrow(PetscBarrier(NULL)); 
            //answer requests
            mySendRequests.resize(0);   mySendRequests.reserve(2*size*size+10); 
            mySendIDs.resize(0);        mySendIDs.reserve(2*size*size+10); 
            myReceiveIDs.resize(0);     myReceiveIDs.reserve(2*size*size+10); 
            std::vector<std::vector<size_t>> DofsToSend;
            std::vector<size_t> RequestVec;
            for (int RequestFrom = 0; RequestFrom < size; RequestFrom++){
                if (RequestNum(RequestFrom, rank)>0){
                    RequestVec.push_back(RequestFrom);
                    mySendRequests.push_back(0); mySendIDs.push_back(rank*size+RequestFrom);
                    DofsToSend.resize(DofsToSend.size()+1); DofsToSend[DofsToSend.size()-1].resize(RequestNum(RequestFrom, rank));
                    getDofForNodesSeries(RequestNodes(RequestFrom,rank), RequestDofs(RequestFrom,rank), DofsToSend[DofsToSend.size()-1]);
                }
            }
            for (size_t i = 0; i < mySendIDs.size(); i++){
                MPI_Isend(DofsToSend[i].data(), DofsToSend[i].size(), my_MPI_SIZE_T, RequestVec[i], mySendIDs[i], PETSC_COMM_WORLD, &mySendRequests[i]);
            }

            //process requests
            std::vector<std::vector<size_t>> DofsReceived(size);
            for (int RequestTo = 0; RequestTo < size; RequestTo++){
                if (RequestNum(rank, RequestTo)>0){
                    mySendRequests.push_back(0); myReceiveIDs.push_back(RequestTo*size+rank);
                    DofsReceived[RequestTo].resize(RequestNum(rank, RequestTo));
                    MPI_Irecv(DofsReceived[RequestTo].data(), DofsReceived[RequestTo].size(), my_MPI_SIZE_T, RequestTo, myReceiveIDs.back(), PETSC_COMM_WORLD, &mySendRequests.back());
                }
            }
            MPI_Waitall(mySendRequests.size(), mySendRequests.data(), myStatus.data());
            PetscCallThrow(PetscBarrier(NULL)); 

            for (int c = 0; c < size; c++){
                if (RequestNum(rank, c)>0){
                    for (size_t idof = 0; idof < DofsReceived[c].size(); idof++){
                        DofNumbering[CurdofStep][RequestDofs(rank,c)[idof]].insert({RequestNodes(rank,c)[idof], DofsReceived[c][idof]});
                        GhostDofs[CurdofStep].push_back(DofsReceived[c][idof]);
                    }
                }
            }
        }

        printStats(CurdofStep);
    }
}

/// @brief Prints statistics regarding degrees of freedom
/// @param curStep Step to print info for
void DofSpace::printStats(size_t curStep){
    size_t GhostTotal = GhostDofs[curStep].size();
    size_t LocalTotal = LocalDofNumRange[curStep][1]-LocalDofNumRange[curStep][0];
    // combine results from all processes
    MPI_Allreduce(MPI_IN_PLACE, &GhostTotal, 1, my_MPI_SIZE_T, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &LocalTotal, 1, my_MPI_SIZE_T, MPI_SUM, PETSC_COMM_WORLD);

    std::string Message = "Degrees of Freedom Added: ";
    for (size_t i = 0; i < nDofs; i++){
        if (DofStep[i]==curStep){
            Message += DofNames[i] + ", ";
        }
    }
    Message += "Total amount of DOFs: "+std::to_string(TotalDofRange[curStep]);
    Logs.PrintSingle(Message+"\n",1);
    Logs.PrintSingle("Total Local dofs:"+std::to_string(LocalTotal)+"\n",2);
    Logs.PrintSingle("Total Ghost dofs:"+std::to_string(GhostTotal)+"\n",2);
    PetscCallThrow(PetscBarrier(NULL)); 
    Message = "Local Dofs: "+std::to_string(LocalDofNumRange[curStep][1]-LocalDofNumRange[curStep][0])+", Ghost Dofs: "+std::to_string(GhostDofs[curStep].size())+"\n";
    Logs.PrintEvery(Message,2);
    PetscCallThrow(PetscBarrier(NULL));
}

void DofSpace::ExportDofSpace(){
    // Create file
    std::string FileName = FilenamePrefix+"/DofSpace.hdf5";
    
    FileAccessProps fapl;
    fapl.add(MPIOFileAccess{PETSC_COMM_WORLD, MPI_INFO_NULL});

    File file(FileName, File::ReadWrite | File::Create | File::Truncate, fapl);

    // Loop through each step and save the data
    for (size_t step = 0; step < maxSteps; step++){
        std::string stepName = "Step_" + std::to_string(step);
        file.createGroup(stepName);
        
        // Save the DofNumbering for this step
        for (size_t dofInd = 0; dofInd < nDofs; dofInd++){
            std::string DofGroupName = stepName + "/" + DofNames[dofInd];
            file.createGroup(DofGroupName);

            std::map<size_t, size_t> m = DofNumbering[step][dofInd];
            size_t dofNums = m.size();
            std::vector<size_t> nodeNumbers;
            std::vector<size_t> dofNumbers;
            
            for(std::map<size_t,size_t>::iterator it = m.begin(); it != m.end(); ++it) {
                nodeNumbers.push_back(it->first);
                dofNumbers.push_back(it->second);
            }

            std::vector<size_t> TotalDofNums(size);
            MPI_Allgather(&dofNums, 1, my_MPI_SIZE_T, TotalDofNums.data(), 1, my_MPI_SIZE_T, PETSC_COMM_WORLD);
            for (size_t i = 0; i < size; i++){
                DataSet dataset = file.createDataSet<size_t>(DofGroupName + "/Core_" + std::to_string(i), DataSpace({TotalDofNums[i],2}));
                if(i == rank){
                    dataset.select({0,0},{dofNums,1} ).write(nodeNumbers);
                    dataset.select({0,1},{dofNums,1} ).write(dofNumbers);
                }
            }
        }
    }

}

/// @brief Obtains degree of freedom indices for the provided node
/// @param Node input: node to obtain the degree of freedom index for (global numbering)
/// @param DofInd input: Index of degree of freedom type to obtain
/// @param Dof_out output: degree of freedom numbering
void DofSpace::getDofForNodes(size_t &Node, size_t DofInd, PetscInt &Dof_out){
    Dof_out = DofNumbering[DofStep[DofInd]][DofInd][Node];
}



/// @brief Obtains degrees of freedom indices for the provided nodes
/// @param Nodes input: nodes to obtain the degree of freedom index for (global numbering)
/// @param DofInd input: Index of degree of freedom type to obtain
/// @param Dofs_out output: degree of freedom numbering (needs to be pre-sized)
void DofSpace::getDofForNodes(std::vector<size_t> &Nodes, size_t DofInd, std::vector<size_t> &Dofs_out){
    for (size_t i = 0; i < Nodes.size(); i++){
        Dofs_out[i] = DofNumbering[DofStep[DofInd]][DofInd][Nodes[i]];
    }
}

/// @brief Obtains degrees of freedom indices for the provided nodes
/// @param Nodes input: nodes to obtain the degree of freedom index for (global numbering)
/// @param DofInd input: Index of degree of freedom type to obtain
/// @param Dofs_out output: degree of freedom numbering (needs to be pre-sized)
void DofSpace::getDofForNodes(std::vector<size_t> &Nodes, size_t DofInd, std::vector<PetscInt> &Dofs_out){
    for (size_t i = 0; i < Nodes.size(); i++){
        Dofs_out[i] = DofNumbering[DofStep[DofInd]][DofInd][Nodes[i]];
    }
}

/// @brief Obtains degrees of freedom indices for the provided nodes
/// @param Nodes input: nodes to obtain the degree of freedom index for (global numbering)
/// @param DofInd input: Index of degree of freedom type to obtain
/// @param Dofs_out output: degree of freedom numbering (needs to be pre-sized)
void DofSpace::getDofForNodes(std::vector<size_t> &Nodes, std::vector<size_t> DofInd, std::vector<PetscInt> &Dofs_out){
    for (size_t j = 0; j < DofInd.size(); j++){
        for (size_t i = 0; i < Nodes.size(); i++){
            Dofs_out[i+j*Nodes.size()] = DofNumbering[DofStep[DofInd[j]]][DofInd[j]][Nodes[i]];
        }
    }
}

/// @brief Obtains degrees of freedom indices for the provided combination of nodes and dof index
/// @param Nodes input: nodes to obtain the degree of freedom index for (global numbering)
/// @param DofInd input: Index of degree of freedom type to obtain
/// @param Dofs_out output: degree of freedom numbering (needs to be pre-sized)
void DofSpace::getDofForNodesSeries(std::vector<size_t> &Nodes, std::vector<size_t> &DofInd, std::vector<size_t> &Dofs_out){
    for (size_t i = 0; i < Nodes.size(); i++){
        Dofs_out[i] = DofNumbering[DofStep[DofInd[i]]][DofInd[i]][Nodes[i]];
    }
}

/// @brief Obtain ghost numbering schemes
/// @param Step input: staggered step
/// @return ghost degree of freedom indices
std::vector<PetscInt> DofSpace::getGhostDofNumbers(size_t Step){
    return GhostDofs[Step];
}