#include "mesh.h"
#include "../InputsOutputs/inputData.h"
#include "../utility/utility.h"

#include <petsc.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <highfive/H5File.hpp>


/// @brief Initialize mesh based on input file
/// @param inputs 
Mesh::Mesh(inputData& inputs) {
    PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
    
    std::string inputName; inputs.GetRequired(inputName,{"mesh","file"});
    inputs.GetRequired(dim,{"mesh","dim"});
    inputs.GetRequired(ipcount1D,{"mesh","ipcount1D"});

    LoadFromFile(inputName);
    Check();
}

/// @brief Initialize mesh based on input file and restart datafile
/// @param inputs 
/// @param data 
Mesh::Mesh(inputData& inputs, SaveDataFile& data){
    //assuming no changes in mesh
    PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));

    std::string inputName; inputs.GetRequired(inputName,{"mesh","file"});
    inputs.GetRequired(dim,{"mesh","dim"});
    inputs.GetRequired(ipcount1D,{"mesh","ipcount1D"});
    LoadFromFile(inputName);
    Check();
}

/// @brief Save parameters required to be safely restarted
/// @param data 
void Mesh::Save(SaveDataFile& data){

}

Mesh::~Mesh() {

}

/// @brief Loads mesh from a provided hdf5 input file
/// @param filename path to the mesh file
void Mesh::LoadFromFile(std::string filename){
    //open mesh input file
    Logs.PrintSingle("Loading mesh from: " + filename + "\n",1);

    HighFive::FileAccessProps fapl;
    fapl.add(HighFive::MPIOFileAccess{PETSC_COMM_WORLD, MPI_INFO_NULL});
    
    HighFive::File file(filename, HighFive::File::ReadOnly, fapl );

    if (size>1){
        HighFive::DataSet CoresForMeshSet = file.getDataSet("/NCores");
        PetscMPIInt NCores; CoresForMeshSet.read(NCores);
        if (NCores != size){
            throw std::invalid_argument("Mesh is not generated for this number of CPU cores");
        }    
    }
    LoadNodesFromFile(file);
    Logs.PrintSingle("Nodes Loaded\n", 2); PetscCallThrow(PetscBarrier(NULL));
    LoadNodeGroupsFromFile(file);
    Logs.PrintSingle("NodeGroups Loaded\n", 2); PetscCallThrow(PetscBarrier(NULL));
    LoadElementGroupsFromFile(file);
    Logs.PrintSingle("Elements Loaded\n", 2); PetscCallThrow(PetscBarrier(NULL));

    //distributing nodes and parrallel ghost vector sync
    Xcoords.AssemblyEnd(); Ycoords.AssemblyEnd();
    if (dim==3){
        Zcoords.AssemblyEnd();
    }
    PetscCallThrow(PetscBarrier(NULL));
    Xcoords.SyncStart(INSERT_VALUES); Ycoords.SyncStart(INSERT_VALUES);
    Xcoords.SyncEnd(INSERT_VALUES); Ycoords.SyncEnd(INSERT_VALUES);
    if (dim==3){
        Zcoords.SyncStart(INSERT_VALUES);
        Zcoords.SyncEnd(INSERT_VALUES);
    }
    PetscCallThrow(PetscBarrier(NULL));
}

void Mesh::LoadNodesFromFile(HighFive::File& file){
    // Get Number of Nodes
    HighFive::DataSet NodeSet = file.getDataSet("/nodes");
    std::vector<size_t> DimNodes = NodeSet.getDimensions();
    size_t NNodes = DimNodes[1];
    Logs.PrintSingle("Number of Nodes: "+std::to_string(NNodes)+"\n",2);

    if (size==1){   //single-core mesh
        Xcoords.SetSizeSequential(NNodes);
        Ycoords.SetSizeSequential(NNodes);
        if (dim==3){
            Zcoords.SetSizeSequential(NNodes);
        }

        std::vector<std::vector<double>> Nodes_stdVec;
        std::vector<PetscInt> Node_Alloc(NNodes); std::iota(Node_Alloc.begin(), Node_Alloc.end(),0);

        NodeSet.read(Nodes_stdVec);
        Xcoords.Set(Node_Alloc, Nodes_stdVec[0], INSERT_VALUES);  Xcoords.AssemblyStart();
        Ycoords.Set(Node_Alloc, Nodes_stdVec[1], INSERT_VALUES);  Ycoords.AssemblyStart();
        if (dim==3){
            Zcoords.Set(Node_Alloc, Nodes_stdVec[2], INSERT_VALUES);  Zcoords.AssemblyStart();
        }

        bool NodesHaveData = file.exist("/node_DataTypes");
        if (NodesHaveData){
            HighFive::DataSet DataNamesSet = file.getDataSet("/node_DataTypes");
            DataNamesSet.read(NodeDataNames);

            HighFive::DataSet NodeDataSet = file.getDataSet("/node_Data");
            std::vector<std::vector<double>> Nodes_DataVec;
            NodeDataSet.read(Nodes_DataVec);

            NodeData.resize(NodeDataNames.size());
            Logs.PrintSingle("Node Data: \n",2);
            for (size_t i = 0; i < NodeDataNames.size(); i++){
                NodeData[i].SetSizeSequential(NNodes);
                NodeData[i].Set(Node_Alloc, Nodes_DataVec[i], INSERT_VALUES); 
                NodeData[i].AssemblyStart(); 
                NodeData[i].AssemblyEnd();
                Logs.PrintSingle("\t"+NodeDataNames[i]+"\n",2);
            }
        }
    } else {
        std::vector<size_t> NodeRange; 
        size_t nNodes; //local nodes
        size_t nGhost; //number of ghost nodes
        std::vector<PetscInt> locGhosts;

        HighFive::DataSet HasGhost = file.getDataSet("/Partition/"+std::to_string(rank)+"/Hasghosts");
        uint8_t hasghost; HasGhost.read(hasghost);

        if (hasghost == true){
            HighFive::DataSet GhostSet = file.getDataSet("/Partition/"+std::to_string(rank)+"/ghosts");
            GhostSet.read(locGhosts);
            nGhost = locGhosts.size();
        } else {
            nGhost = 0;
            locGhosts.resize(0);
        }

        HighFive::DataSet NodeRangeSet = file.getDataSet("/Partition/"+std::to_string(rank)+"/noderange");
        NodeRangeSet.read(NodeRange);
        nNodes = NodeRange[1] - NodeRange[0]+1;

        Xcoords.SetSize(nNodes, NNodes, locGhosts);
        Ycoords.SetSize(nNodes, NNodes, locGhosts);
        if (dim==3){
            Zcoords.SetSize(nNodes, NNodes, locGhosts);
        }

        std::vector<std::vector<double>> Nodes_stdVec;
        std::vector<PetscInt> Node_Alloc(nNodes); std::iota(Node_Alloc.begin(), Node_Alloc.end(),NodeRange[0]);
        std::vector<size_t> Node_Alloc2(nNodes); std::iota(Node_Alloc2.begin(), Node_Alloc2.end(),NodeRange[0]);
        NodeSet.select(Node_Alloc2).read(Nodes_stdVec);
        
        Logs.PrintEvery("Owned: "+std::to_string(nNodes)+", ghosts: "+std::to_string(nGhost)+"\n",3);
        Xcoords.Set(Node_Alloc, Nodes_stdVec[0], INSERT_VALUES);  Xcoords.AssemblyStart();
        Ycoords.Set(Node_Alloc, Nodes_stdVec[1], INSERT_VALUES);  Ycoords.AssemblyStart();
        if (dim==3){
            Zcoords.Set(Node_Alloc, Nodes_stdVec[2], INSERT_VALUES);  Zcoords.AssemblyStart();
        }

        bool NodesHaveData = file.exist("/node_DataTypes");
        if (NodesHaveData){
            HighFive::DataSet DataNamesSet = file.getDataSet("/node_DataTypes");
            DataNamesSet.read(NodeDataNames);

            HighFive::DataSet NodeDataSet = file.getDataSet("/node_Data");
            std::vector<std::vector<double>> Nodes_DataVec;
            NodeDataSet.select(Node_Alloc2).read(Nodes_DataVec);

            NodeData.resize(NodeDataNames.size());
            Logs.PrintSingle("Node Data: \n",2);
            for (size_t i = 0; i < NodeDataNames.size(); i++){
                NodeData[i].SetSize(nNodes, NNodes, locGhosts);
                NodeData[i].Set(Node_Alloc, Nodes_DataVec[i], INSERT_VALUES); 
                NodeData[i].AssemblyStart(); 
                NodeData[i].AssemblyEnd();
                NodeData[i].SyncStart(INSERT_VALUES);
                NodeData[i].SyncEnd(INSERT_VALUES);
                Logs.PrintSingle("\t"+NodeDataNames[i]+"\n",2);
            }
        }
    }
}

void Mesh::LoadNodeGroupsFromFile(HighFive::File& file){
    HighFive::DataSet nodegroupnamesSet = file.getDataSet("/nodegroupnames");
    nodegroupnamesSet.read(NodeGroupNames);
    size_t NNGroups = NodeGroupNames.size();

    if (size==1){   //single-core mesh
        Logs.PrintSingle("NodeGroups:\n",2);
        NodeGroups.resize(NNGroups);
        for (size_t i = 0; i < NNGroups; i++){
            HighFive::DataSet CurrentGroupSet = file.getDataSet("/nodegroups/"+NodeGroupNames[i]);
            std::vector<size_t> NodesInGroup;
            CurrentGroupSet.read(NodesInGroup);
            NodeGroups[i].AddNodes(NodesInGroup);
            NodeGroups[i].SetName(NodeGroupNames[i]);
            Logs.PrintSingle("\t"+NodeGroups[i].Name +": "+std::to_string(NodeGroups[i].NNodes)+"\n",2);
        }
    } else {
        Logs.PrintSingle("NodeGroups:\n",2);
        NodeGroups.resize(NNGroups);
        HighFive::DataSet HasNodesGroupSet = file.getDataSet("/Partition/"+std::to_string(rank)+"/hasnodegroup");
        std::vector<uint8_t> HasNodesGroup; HasNodesGroupSet.read(HasNodesGroup);

        HighFive::DataSet CurrentGroupSet = file.getDataSet("/Partition/"+std::to_string(rank)+"/nodegrouprange");
        std::vector<std::vector<size_t>> NodeGroupRange; CurrentGroupSet.read(NodeGroupRange);
        for (size_t i = 0; i < NNGroups; i++){
            NodeGroups[i].SetName(NodeGroupNames[i]);
            HighFive::DataSet CurrentGroupSet = file.getDataSet("/nodegroups/"+NodeGroupNames[i]);
            if (HasNodesGroup[i]==true){
                std::vector<size_t> NodeRangeVec(NodeGroupRange[1][i]-NodeGroupRange[0][i]+1); std::iota(NodeRangeVec.begin(), NodeRangeVec.end(),NodeGroupRange[0][i]);
                std::vector<size_t> NodesInGroup;
                CurrentGroupSet.select(NodeRangeVec).read(NodesInGroup);    
                NodeGroups[i].AddNodes(NodesInGroup);            
            }
            std::vector<size_t> GroupSize = CurrentGroupSet.getDimensions();
            Logs.PrintSingle("\t"+NodeGroups[i].Name +": "+std::to_string(GroupSize[0])+"\n",2);
        }
    }
}

void Mesh::LoadElementGroupsFromFile(HighFive::File& file){
    HighFive::DataSet elementgroupnamesSet = file.getDataSet("/elementgroupnames");
    elementgroupnamesSet.read(ElemGroupNames);
    size_t NEGroups = ElemGroupNames.size();

    HighFive::DataSet elementgrouptypeSet = file.getDataSet("/elementgrouptypes");
    std::vector<std::string> ElemGroupTypes;
    elementgrouptypeSet.read(ElemGroupTypes);

    if (size==1){   //single-core mesh
        Logs.PrintSingle("ElementGroups:\n",2);
        ElementGroups.resize(NEGroups);
        for (size_t i = 0; i < NEGroups; i++){
            HighFive::DataSet CurrentGroupSet = file.getDataSet("/elementgroups/"+ElemGroupNames[i]);
            std::vector<std::vector<size_t>> ElemInGroup;
            ElementGroups[i].SetName(ElemGroupNames[i]);
            ElementGroups[i].setType(ElemGroupTypes[i], dim, ipcount1D);
            CurrentGroupSet.read(ElemInGroup);

            if (ElementGroups[i].BaseElem->requiresData){
                HighFive::DataSet CurrentGroupDataSet = file.getDataSet("/elementgroups/"+ElemGroupNames[i]+"_Data");
                std::vector<std::vector<std::vector<double>>> ElemData;
                CurrentGroupDataSet.read(ElemData);
                ElementGroups[i].AddElems(ElemInGroup, ElemData);
            } else {
                ElementGroups[i].AddElems(ElemInGroup);
            }

            Logs.PrintSingle("\t"+ElementGroups[i].Name +": "+std::to_string(ElementGroups[i].NElems)+"  ("+ElementGroups[i].BaseElem->Name+")\n",2);
        }
    } else { // parrallel mesh
        Logs.PrintSingle("ElementGroups:\n",2); PetscCallThrow(PetscBarrier(NULL));

        ElementGroups.resize(NEGroups);
        HighFive::DataSet HasNElementGroupSet = file.getDataSet("/Partition/"+std::to_string(rank)+"/Haselems");
        std::vector<uint8_t> HasElemGroup; HasNElementGroupSet.read(HasElemGroup);

        HighFive::DataSet ElemGroupSet = file.getDataSet("/Partition/"+std::to_string(rank)+"/elemrange");
        std::vector<std::vector<size_t>> ElemGroupRange; ElemGroupSet.read(ElemGroupRange);

        for (size_t i = 0; i < NEGroups; i++){
            ElementGroups[i].SetName(ElemGroupNames[i]);

            Logs.PrintSingle("\t"+ElemGroupNames[i]+"\n",2); PetscCallThrow(PetscBarrier(NULL));

            ElementGroups[i].setType(ElemGroupTypes[i], dim, ipcount1D);

            //Logs.PrintSingle("Types Set\n",2); PetscCallThrow(PetscBarrier(NULL));
            if (HasElemGroup[i]==true){
                std::vector<size_t> ElemRangeVec(ElemGroupRange[1][i]-ElemGroupRange[0][i]+1); std::iota(ElemRangeVec.begin(), ElemRangeVec.end(),ElemGroupRange[0][i]);
                
                HighFive::DataSet CurrentGroupSet = file.getDataSet("/elementgroups/"+ElemGroupNames[i]);
                std::vector<std::vector<size_t>> ElemsInGroup;
                CurrentGroupSet.select(ElemRangeVec).read(ElemsInGroup);    
                
                if (ElementGroups[i].BaseElem->requiresData){
                    HighFive::DataSet CurrentGroupDataSet = file.getDataSet("/elementgroups/"+ElemGroupNames[i]+"_Data");
                    std::vector<std::vector<std::vector<double>>> ElemData;
                    CurrentGroupDataSet.select(ElemRangeVec).read(ElemData);
                    ElementGroups[i].AddElems(ElemsInGroup, ElemData);
                } else {
                    ElementGroups[i].AddElems(ElemsInGroup);
                }     
            }
        }
    }
}

/// @brief Calculate the area of every element group, verifying all elements can be evaluated (and allowing for manually confirming the area is correct)
void Mesh::Check(){
    Logs.PrintSingle("Checking mesh\n",1);
    size_t NGroups = ElementGroups.size();
    Eigen::VectorXd LocalArea(NGroups), AllArea(NGroups);

    std::vector<size_t> Elem_Nodes;
    Eigen::MatrixXd coordsNodes;
    std::vector<Eigen::RowVectorXd> N;
    std::vector<Eigen::MatrixXd> G;
    std::vector<double> w;

    //Loop over all element groups
    LocalArea.setZero();
    for (size_t grpInd = 0; grpInd < NGroups; grpInd++){
        size_t NperElem = ElementGroups[grpInd].NNodes_per_elem;
        Elem_Nodes.resize(NperElem);
        coordsNodes.resize(NperElem, dim);
        N.resize(ElementGroups[grpInd].BaseElem->ipcount); for (size_t i = 0; i < N.size(); i++) N[i].resize(NperElem);
        G.resize(ElementGroups[grpInd].BaseElem->ipcount); for (size_t i = 0; i < G.size(); i++) G[i].resize(ElementGroups[grpInd].BaseElem->dim, NperElem);
        w.resize(ElementGroups[grpInd].BaseElem->ipcount);

        for (size_t el = 0; el < ElementGroups[grpInd].NElems; el++){
            //get element data, and perform numerical integration (summing weights)
            getShapeGrads(grpInd, el, Elem_Nodes, N, G, w);
            for (double wi : w) LocalArea(grpInd) += wi;
        }
        MPI_Allreduce(&LocalArea[grpInd], &AllArea[grpInd], 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
        PetscCallThrow(PetscBarrier(NULL)); 

        //Print resulting areas
        std::stringstream ss; ss << "\t" << ElemGroupNames[grpInd] << ": " << AllArea(grpInd) << "  (Local: " << LocalArea(grpInd) << " )\n";
        Logs.PrintSingle(ss.str(),1);
    }
}

/// @brief Provide the shape functions, shape function gradients, and integration weights. output vectors need to be pre-sized
/// @param egroup Element group index to evaluate for
/// @param nEl Element index (local numbering) to provide outputs for
/// @param Nodes output: Node indices contained within this element
/// @param Nout  output: Element shape functions within integration points, accessed through N[ip] for shape function rowvector of size nNodes
/// @param Gout output: Gradient of shape functions, accessed through G[ip] for gradient matrix of size (dim,nNodes)
/// @param wout output: Integration weights
void Mesh::getShapeGrads(size_t egroup, size_t nEl, std::vector<size_t>& Nodes, std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout){
    Eigen::MatrixXd CoordsNodes(ElementGroups[egroup].NNodes_per_elem,dim);    
    
    GetNodesForElem(Nodes, egroup, nEl);
    GetCoordsForNodes(CoordsNodes, Nodes);
    
    if (ElementGroups[egroup].BaseElem->requiresData && ElementGroups[egroup].BaseElem->requiresNodeData){
        Eigen::VectorXd NodeDataVec(ElementGroups[egroup].NNodes_per_elem);
        NodeData[0].GetValues(Nodes, NodeDataVec);
        ElementGroups[egroup].BaseElem->getShapeGrads(Nout, Gout, wout, CoordsNodes, ElementGroups[egroup].ElemData[nEl], NodeDataVec);
    } else if (ElementGroups[egroup].BaseElem->requiresData){
        ElementGroups[egroup].BaseElem->getShapeGrads(Nout, Gout, wout, CoordsNodes, ElementGroups[egroup].ElemData[nEl]);
    } else {
        ElementGroups[egroup].BaseElem->getShapeGrads(Nout, Gout, wout, CoordsNodes);
    }

}

/// @brief Provide the shape functions, 1st and 2nd shape function gradients, and integration weights. output vectors need to be pre-sized
/// @param egroup Element group index to evaluate for
/// @param nEl Element index (local numbering) to provide outputs for
/// @param Nodes output: Node indices contained within this element
/// @param Nout  output: Element shape functions within integration points, accessed through N[ip] for shape function rowvector of size nNodes
/// @param Gout output: Gradient of shape functions, accessed through G[ip] for gradient matrix of size (dim,nNodes)
/// @param G2out output: Second gradient of shape functions, accessed through G2[ip] for gradient matrix of size (xx/yy/xy,nNodes)
/// @param wout output: Integration weights
void Mesh::getShapeGrads(size_t egroup, size_t nEl, std::vector<size_t>& Nodes, std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout){
    Eigen::MatrixXd CoordsNodes(ElementGroups[egroup].NNodes_per_elem,dim);    
    
    GetNodesForElem(Nodes, egroup, nEl);
    GetCoordsForNodes(CoordsNodes, Nodes);
    
    if (ElementGroups[egroup].BaseElem->requiresData && ElementGroups[egroup].BaseElem->requiresNodeData){
        Eigen::VectorXd NodeDataVec(ElementGroups[egroup].NNodes_per_elem);
        NodeData[0].GetValues(Nodes, NodeDataVec);
        ElementGroups[egroup].BaseElem->getShapeGrads2(Nout, Gout, G2out, wout, CoordsNodes, ElementGroups[egroup].ElemData[nEl], NodeDataVec);
    } else if (ElementGroups[egroup].BaseElem->requiresData){
        ElementGroups[egroup].BaseElem->getShapeGrads2(Nout, Gout, G2out, wout, CoordsNodes, ElementGroups[egroup].ElemData[nEl]);
    } else {
        ElementGroups[egroup].BaseElem->getShapeGrads2(Nout, Gout, G2out, wout, CoordsNodes);
    }

}

/// @brief Provides the coordinates associated with each integration point
/// @param egroup Element group index to evaluate for
/// @param nEl Element index (local numbering) to provide coordinates for
/// @param coordsIP output: Matrix containing integration point coordinates, accessed as coordsIP(ip,dim))
/// @param coordsNodes input: Nodal coordinates (NOT RECALCULATED)
void Mesh::getCoordsIP(size_t egroup, size_t nEl, Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes){
    if (ElementGroups[egroup].BaseElem->requiresData && ElementGroups[egroup].BaseElem->requiresNodeData){
        std::vector<size_t> Nodes; Nodes.resize(ElementGroups[egroup].NNodes_per_elem);
        GetNodesForElem(Nodes, egroup, nEl);
        Eigen::VectorXd NodeDataVec(ElementGroups[egroup].NNodes_per_elem);
        NodeData[0].GetValues(Nodes, NodeDataVec);
        ElementGroups[egroup].BaseElem->getCoordsIP(coordsIP, coordsNodes, ElementGroups[egroup].ElemData[nEl], NodeDataVec);
    } else if (ElementGroups[egroup].BaseElem->requiresData){
        ElementGroups[egroup].BaseElem->getCoordsIP(coordsIP, coordsNodes, ElementGroups[egroup].ElemData[nEl]);
    } else {
        ElementGroups[egroup].BaseElem->getCoordsIP(coordsIP, coordsNodes);
    }   
}

void Mesh::getNormals(size_t egroup, size_t nEl, std::vector<Eigen::VectorXd> &normals){
    std::vector<size_t> Nodes(ElementGroups[egroup].NNodes_per_elem);
    Eigen::MatrixXd CoordsNodes(ElementGroups[egroup].NNodes_per_elem,dim);    
    
    GetNodesForElem(Nodes, egroup, nEl);
    GetCoordsForNodes(CoordsNodes, Nodes);

    if (ElementGroups[egroup].BaseElem->requiresData && ElementGroups[egroup].BaseElem->requiresNodeData){
        std::vector<size_t> Nodes; Nodes.resize(ElementGroups[egroup].NNodes_per_elem);
        GetNodesForElem(Nodes, egroup, nEl);
        Eigen::VectorXd NodeDataVec(ElementGroups[egroup].NNodes_per_elem);
        NodeData[0].GetValues(Nodes, NodeDataVec);
        ElementGroups[egroup].BaseElem->getNormals(normals, CoordsNodes, ElementGroups[egroup].ElemData[nEl], NodeDataVec);
    } else if (ElementGroups[egroup].BaseElem->requiresData){
        ElementGroups[egroup].BaseElem->getNormals(normals, CoordsNodes, ElementGroups[egroup].ElemData[nEl]);
    } else {
        ElementGroups[egroup].BaseElem->getNormals(normals, CoordsNodes);
    }   
}

/// @brief Get shape functions at corners of element for exporting nodal quantities
/// @param egroup Element group index to evaluate for
/// @param nEl Element index (local numbering) to provide coordinates for
/// @param Nodes UNUSED
/// @param Nout Export shape functions, array accessed as N[point] to produce a rowvector of size nNodes
void Mesh::getExportShape(size_t egroup, size_t nEl, std::vector<size_t>& Nodes, std::vector<Eigen::RowVectorXd> &Nout){
    if (ElementGroups[egroup].BaseElem->requiresData && ElementGroups[egroup].BaseElem->requiresNodeData){
        Eigen::VectorXd NodeDataVec(ElementGroups[egroup].NNodes_per_elem);
        NodeData[0].GetValues(Nodes, NodeDataVec);
        ElementGroups[egroup].BaseElem->getExportShape(Nout, ElementGroups[egroup].ElemData[nEl], NodeDataVec);
    } else if (ElementGroups[egroup].BaseElem->requiresData){
        ElementGroups[egroup].BaseElem->getExportShape(Nout, ElementGroups[egroup].ElemData[nEl]);
    } else {
        ElementGroups[egroup].BaseElem->getExportShape(Nout);
    }   
}

/// @brief Obtain nodes associated to the current element
/// @param outNodes output: Nodes for current element
/// @param NGroup Element group index to evaluate for
/// @param NElem Element index (local numbering) to evaluate for
void Mesh::GetNodesForElem(std::vector<size_t>& outNodes, size_t NGroup, size_t NElem){
    for (size_t i = 0; i < ElementGroups[NGroup].NNodes_per_elem; i++){
        outNodes[i] = ElementGroups[NGroup].Elems[NElem][i];
    }
}

/// @brief Obtain the coordinates for the provided node
/// @param outCoordsNodes output: Coordinates, given as Coords(node,dim)
/// @param Nodes input: Node indices
void Mesh::GetCoordsForNodes(Eigen::MatrixXd& outCoordsNodes, std::vector<size_t>& Nodes){
    Eigen::VectorXd X, Y, Z; X.resize(Nodes.size()); Y.resize(Nodes.size()); Z.resize(Nodes.size());

    Xcoords.GetValues(Nodes, X); outCoordsNodes.col(0) = X;
    Ycoords.GetValues(Nodes, Y); outCoordsNodes.col(1) = Y;
    if (dim==3){
        Zcoords.GetValues(Nodes, Z); outCoordsNodes.col(2) = Z;
    }
}

void Mesh::GetCoordsForNodes(Eigen::VectorXd& outCoordsNode, size_t& Node){
    double x,y,z;
    x = Xcoords.GetValue(Node); outCoordsNode(0) = x;
    y = Ycoords.GetValue(Node); outCoordsNode(1) = y;
    if (dim==3){
        z = Zcoords.GetValue(Node); outCoordsNode(2) = z;
    } 
}

/// @brief Obtain the coordinates for the provided node
/// @param XCoordsOut output: Vector with x coordinates
/// @param YCoordsOut output: vector with y coordinates
/// @param Nodes input: Node indices
void Mesh::GetCoordsForNodes(std::vector<double>& XCoordsOut, std::vector<double>& YCoordsOut, std::vector<size_t>& Nodes){
    Xcoords.GetValues(Nodes, XCoordsOut);
    Ycoords.GetValues(Nodes, YCoordsOut);
    if (dim==3){
        std::invalid_argument("Mesh GetCoordsForNodes 2 function is called with 2 arguments for a 3d mesh");
    }
}

/// @brief Obtain the coordinates for the provided node
/// @param XCoordsOut output: Vector with x coordinates
/// @param YCoordsOut output: vector with y coordinates
/// @param Nodes input: Node indices
void Mesh::GetCoordsForNodes(Eigen::VectorXd& XCoordsOut, Eigen::VectorXd& YCoordsOut, std::vector<size_t>& Nodes){
    Xcoords.GetValues(Nodes, XCoordsOut);
    Ycoords.GetValues(Nodes, YCoordsOut);
    if (dim==3){
        std::invalid_argument("Mesh GetCoordsForNodes 3 function is called with 2 arguments for a 3d mesh");
    }
}

/// @brief Obtain the coordinates for the provided node
/// @param XCoordsOut output: Vector with x coordinates
/// @param YCoordsOut output: vector with y coordinates
/// @param ZCoordsOut output: vector with z coordinates
/// @param Nodes input: Node indices
void Mesh::GetCoordsForNodes(std::vector<double>& XCoordsOut, std::vector<double>& YCoordsOut, std::vector<double>& ZCoordsOut, std::vector<size_t>& Nodes){
    if (dim==2){
        std::invalid_argument("Mesh GetCoordsForNodes 2 function is called with 2 arguments for a 3d mesh");
    }
    Xcoords.GetValues(Nodes, XCoordsOut);
    Ycoords.GetValues(Nodes, YCoordsOut);
    Zcoords.GetValues(Nodes, ZCoordsOut);
}

/// @brief Obtain the coordinates for the provided node
/// @param XCoordsOut output: Vector with x coordinates
/// @param YCoordsOut output: vector with y coordinates
/// @param ZCoordsOut output: vector with z coordinates
/// @param Nodes input: Node indices
void Mesh::GetCoordsForNodes(Eigen::VectorXd& XCoordsOut, Eigen::VectorXd& YCoordsOut, Eigen::VectorXd& ZCoordsOut, std::vector<size_t>& Nodes){
    if (dim==2){
        std::invalid_argument("Mesh GetCoordsForNodes 3 function is called with 2 arguments for a 3d mesh");
    }
    Xcoords.GetValues(Nodes, XCoordsOut);
    Ycoords.GetValues(Nodes, YCoordsOut);
    Zcoords.GetValues(Nodes, ZCoordsOut);
}

/// @brief Converts an element group name to an associated element group index
/// @param GroupName Element group name
/// @return Element group index
size_t Mesh::GetElementGroupIdx(std::string GroupName){
    auto it = std::find(ElemGroupNames.begin(), ElemGroupNames.end(), GroupName);
    size_t index;
    if (it != ElemGroupNames.end()) {
        index = it - ElemGroupNames.begin();
    } else {
        throw std::invalid_argument("Element group \"" + GroupName + "\" is not defined");
    }
    return index;
}

/// @brief Converts a node group name to an associated node group index
/// @param GroupName Node group name
/// @return Node group index
size_t Mesh::GetNodeGroupIdx(std::string GroupName){
    auto it = std::find(NodeGroupNames.begin(), NodeGroupNames.end(), GroupName);
    size_t index;
    if (it != NodeGroupNames.end()) {
        index = it - NodeGroupNames.begin();
    } else {
        throw std::invalid_argument("Node group \"" + GroupName + "\" is not defined");
    }
    return index;
}

size_t Mesh::GetExportMeshCoords(size_t eg, std::vector<std::vector<double>>& XCoords, std::vector<std::vector<double>>& YCoords, std::vector<std::vector<double>>& ZCoords, std::vector<std::vector<double>>& XIPCoords, std::vector<std::vector<double>>& YIPCoords, std::vector<std::vector<double>>& ZIPCoords){
    //get nodal and integration-point coordinates
    size_t nNodes = ElementGroups[eg].NNodes_per_elem;
    size_t nip = ElementGroups[eg].BaseElem->ipcount;
    XCoords.resize(ElementGroups[eg].NElems);
    YCoords.resize(ElementGroups[eg].NElems);
    XIPCoords.resize(ElementGroups[eg].NElems);
    YIPCoords.resize(ElementGroups[eg].NElems);
    if (dim==3){
        ZCoords.resize(ElementGroups[eg].NElems);
        ZIPCoords.resize(ElementGroups[eg].NElems);      
    }

    std::vector<Eigen::RowVectorXd> NExport(ElementGroups[eg].BaseElem->NExport.size()); 
    for (size_t i = 0; i < NExport.size(); i++) NExport[i].resize(nNodes);

    for (size_t i = 0; i < ElementGroups[eg].NElems; i++){
        std::vector<size_t> el(nNodes);
        Eigen::VectorXd XNodes(nNodes), YNodes(nNodes), ZNodes(nNodes);

        GetNodesForElem(el, eg, i);
        getExportShape(eg, i, el, NExport);
        Eigen::MatrixXd coordsNodes(nNodes, dim), coordsIP(nip, dim);
        GetCoordsForNodes(coordsNodes, el);
        if (dim==2){
            GetCoordsForNodes(XNodes, YNodes, el);
        } else if (dim==3){
           GetCoordsForNodes(XNodes, YNodes, ZNodes, el);
        }
        getCoordsIP(eg, i, coordsIP, coordsNodes);        

        XCoords[i].resize(NExport.size());
        YCoords[i].resize(NExport.size());
        XIPCoords[i].resize(nip);
        YIPCoords[i].resize(nip); 

        for (size_t j = 0; j < NExport.size(); j++){
            XCoords[i][j] = NExport[j]*XNodes;
            YCoords[i][j] = NExport[j]*YNodes;
        }
        
        for (size_t j = 0; j < nip; j++){
            XIPCoords[i][j] = coordsIP(j,0);
            YIPCoords[i][j] = coordsIP(j,1);
        }

        if (dim==3){
            ZCoords[i].resize(NExport.size());
            ZIPCoords[i].resize(nip); 

            for (size_t j = 0; j < NExport.size(); j++){
                ZCoords[i][j] = NExport[j]*ZNodes;
            }
            
            for (size_t j = 0; j < nip; j++){
                ZIPCoords[i][j] = coordsIP(j,2);
            }       
        }
    }
    return NExport.size();
}
