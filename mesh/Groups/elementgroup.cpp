#include "elementgroup.h"
#include <iostream>

/// @brief Initialization, sets group as being completely empty
elementgroup::elementgroup(){
    Elems.resize(0);
    NElems = 0;
    Name = "";
};

elementgroup::~elementgroup(){
    delete BaseElem;
};

/// @brief Adds new elements to the element group
/// @param ElemsToAdd Elements to add [Element index][Node index]
void elementgroup::AddElems(std::vector<std::vector<size_t>> ElemsToAdd){
    size_t nElem;
    if (ElemsToAdd.size()==0) {
        nElem = 0;
    } else {
        nElem = ElemsToAdd[0].size();
    }

    size_t NElems_Old = NElems;
    NElems += nElem;
    
    Elems.resize(NElems);
    for (size_t i = 0; i < nElem; i++){
        Elems[NElems_Old+i].resize(NNodes_per_elem);
        for (size_t j = 0; j < NNodes_per_elem; j++){
            Elems[NElems_Old+i][j] = ElemsToAdd[j][i];
        }
    }
 };

 /// @brief Adds new elements with elemental data to the element group
 /// @param ElemsToAdd Elements to add [Element index][Node index]
 /// @param ElemDataToAdd Element data [element index][data_dim_1][data_dim_2]
 void elementgroup::AddElems(std::vector<std::vector<size_t>> ElemsToAdd, std::vector<std::vector<std::vector<double>>> ElemDataToAdd){
    size_t nElem;
    if (ElemsToAdd.size()==0) {
        nElem = 0;
    } else {
        nElem = ElemsToAdd[0].size();
    }

    size_t NElems_Old = NElems;
    NElems += nElem;
    
    Elems.resize(NElems);
    ElemData.resize(NElems);
    size_t DataSize = ElemDataToAdd[0].size();
    for (size_t i = 0; i < nElem; i++){
        Elems[NElems_Old+i].resize(NNodes_per_elem);
        for (size_t j = 0; j < NNodes_per_elem; j++){
            Elems[NElems_Old+i][j] = ElemsToAdd[j][i];
        }
        ElemData[NElems_Old+i].resize(DataSize,DataSize);
        for (size_t j = 0; j < DataSize; j++){
            for (size_t k = 0; k < DataSize; k++){
                ElemData[NElems_Old+i](j,k) = ElemDataToAdd[k][j][i];
            }
        }
    }

 }

/// @brief Sets name of element group
/// @param Name_New New name to use
void elementgroup::SetName(std::string Name_New){
    Name = Name_New;
 }

/// @brief Sets and initializes the parametric element type used for this group
/// @param Type_Name Identifier for the parametric element type to use
/// @param dim UNUSED
/// @param ipcount1D number of integration points per dimension
void elementgroup::setType(std::string Type_Name, int dim, int ipcount1D){
    BaseElem = CreateElem(Type_Name, ipcount1D);
    NNodes_per_elem = BaseElem->NodeCount;
}

/// @brief Get all nodes (global numbering) associated with the elements in this group (duplicates are removed)
/// @return Vector containing all unique nodes
std::vector<size_t> elementgroup::GetUniqueNodes(){
    std::vector<size_t> AllNodes(NNodes_per_elem*NElems);
    for (size_t i = 0; i < NElems; i++){
        for (size_t j = 0; j < NNodes_per_elem; j++){
            AllNodes[i*NNodes_per_elem+j] = Elems[i][j];
        } 
    }
    std::sort(AllNodes.begin(), AllNodes.end());
    auto last = std::unique(AllNodes.begin(), AllNodes.end());
    AllNodes.erase(last, AllNodes.end());
    return AllNodes;
}