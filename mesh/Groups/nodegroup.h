#ifndef NODEGROUP_H
#define NODEGROUP_H

#include <vector>
#include <Eigen/Dense>


/// @brief Group containing indices of nodes associated with set regions
class nodegroup {
    public:
        nodegroup();
        ~nodegroup();
        void AddNodes(std::vector<size_t> ToAdd);
        void SetName(std::string Name);
        
        std::vector<size_t> Nodes; //Node indices (global numdering) associated with this group
        size_t NNodes;  //total number of odes within this group
        std::string Name; //Name of this group (indicates its region)
    private:

};

#endif