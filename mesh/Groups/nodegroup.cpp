#include "nodegroup.h"

nodegroup::nodegroup(){
    Nodes.resize(0);
    NNodes = 0;
    Name = "";
};

nodegroup::~nodegroup(){

 }

 /// @brief Adds nodes to the nodegroup and checks for duplicates
 /// @param ToAdd Nodes (in global numbering) to be added
 void nodegroup::AddNodes(std::vector<size_t> ToAdd){
    size_t NNodes_Old = NNodes;
    NNodes = NNodes + ToAdd.size();
    Nodes.resize(NNodes);
    for (size_t i = 0; i < ToAdd.size(); i++){
        Nodes[NNodes_Old+i] = ToAdd[i];
    }

 }

 /// @brief Sets the name of this node group
 /// @param Name_New Name to be used for group
 void nodegroup::SetName(std::string Name_New){
    Name = Name_New;
 }