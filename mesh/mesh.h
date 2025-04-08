#ifndef MESH_H
#define MESH_H

#include <vector>
#include <Eigen/Dense>
#include <petscvec.h>

#include "../utility/utility.h"
#include "../InputsOutputs/inputData.h"
#include "Groups/nodegroup.h"
#include "Groups/elementgroup.h"
#include "../InputsOutputs/SaveData.h"

/// @brief Mesh used for finite element simulation, holding nodal coordinates, and arrays for mesh and element groups. Also provides shape functions for elements
class Mesh {
    public:
        Mesh(inputData& inputs);
        Mesh(inputData& inputs, SaveDataFile& data);
        void Save(SaveDataFile& data);
        ~Mesh();

        void LoadFromFile(std::string filename);
        void Check();

        void GetNodesForElem(std::vector<size_t>& outNodes, size_t NGroup, size_t NElem);
        void getShapeGrads(size_t egroup, size_t nEl, std::vector<size_t>& Nodes, std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout);
        void getShapeGrads(size_t egroup, size_t nEl, std::vector<size_t>& Nodes, std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout);
        void getExportShape(size_t egroup, size_t nEl, std::vector<size_t>& Nodes, std::vector<Eigen::RowVectorXd> &Nout);
        void getNormals(size_t egroup, size_t nEl, std::vector<Eigen::VectorXd> &normals);
        void GetCoordsForNodes(Eigen::MatrixXd& outCoordsNodes, std::vector<size_t>& Nodes);
        void GetCoordsForNodes(Eigen::VectorXd& outCoordsNode, size_t& Node);
        void GetCoordsForNodes(std::vector<double>& XCoordsOut, std::vector<double>& YCoordsOut, std::vector<size_t>& Nodes);
        void GetCoordsForNodes(Eigen::VectorXd& XCoordsOut, Eigen::VectorXd& YCoordsOut, std::vector<size_t>& Nodes);
        void GetCoordsForNodes(std::vector<double>& XCoordsOut, std::vector<double>& YCoordsOut, std::vector<double>& ZCoordsOut, std::vector<size_t>& Nodes);
        void GetCoordsForNodes(Eigen::VectorXd& XCoordsOut, Eigen::VectorXd& YCoordsOut, Eigen::VectorXd& ZCoordsOut, std::vector<size_t>& Nodes);
        void getCoordsIP(size_t egroup, size_t nEl, Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes);
        size_t GetElementGroupIdx(std::string GroupName);
        size_t GetNodeGroupIdx(std::string GroupName);

        size_t GetExportMeshCoords(size_t egroup, std::vector<std::vector<double>>& XCoords, std::vector<std::vector<double>>& YCoords, std::vector<std::vector<double>>& ZCoords, std::vector<std::vector<double>>& XIPCoords, std::vector<std::vector<double>>& YIPCoords, std::vector<std::vector<double>>& ZIPCoords);

        uint     dim; //mesh dimension
        uint     ipcount1D; //number of integration points per dimension
        std::vector<std::string> ElemGroupNames; //array containing all the element group names
        std::vector<std::string> NodeGroupNames; //array containing all the node group names

        GhostVector Xcoords; //Nodal X coordinates
        GhostVector Ycoords; //Nodal Y coordinates
        GhostVector Zcoords; //Nodal Z coordinates (only used when dim==3)

        std::vector<std::string> NodeDataNames;
        std::vector<GhostVector> NodeData;

        Eigen::RowVector2i MyNodeRange; //lower and upper bounds of node indices contained within this core

        std::vector<nodegroup> NodeGroups; //array containing all the element groups
        std::vector<elementgroup> ElementGroups;//array containing all the node groups
    private:
        PetscMPIInt    rank, size;

        void LoadNodesFromFile(HighFive::File& file);
        void LoadNodeGroupsFromFile(HighFive::File& file);
        void LoadElementGroupsFromFile(HighFive::File& file);
};

#endif