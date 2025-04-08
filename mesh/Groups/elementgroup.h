#ifndef ELEMENTGROUP_H
#define ELEMENTGROUP_H

#include "../ElementTypes/TypeRegister.h"

#include <vector>
#include <Eigen/Dense>

/// @brief Group containing elements (element index in local numbering, nodes pertaining to them in global numbering)
class elementgroup {
    public:
        elementgroup();
        ~elementgroup();
        void AddElems(std::vector<std::vector<size_t>> ElemsToAdd);
        void AddElems(std::vector<std::vector<size_t>> ElemsToAdd, std::vector<std::vector<std::vector<double>>> ElemData);
        void SetName(std::string Name);
        void setType(std::string Type_Name, int dim, int ipcount1D);
        std::vector<size_t> GetUniqueNodes();

        std::vector<std::vector<size_t>> Elems; // Elements [elemIndex][Elem_Nodes] with elemIndex in local numbering, and node indices in global
        std::vector<Eigen::MatrixXd> ElemData;  // Element-specific data required for evaluating shape functions
        size_t NElems;  //total number of elements contained within this group for this specific core
        size_t NNodes_per_elem; //number of nodes per element (required to be constant throughout the element group)
        std::string Name; //Identifier of the element group
        BaseElemType* BaseElem; //Parametric space element used to evaluate shape functions
        
    private:

};

#endif