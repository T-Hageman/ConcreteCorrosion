#ifndef EXTERNALFORCEMODEL_H
#define EXTERNALFORCEMODEL_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"

/// @brief Allows for external forcing to be applied, applying an external force [F_1 F_2 F_3] to degrees of freedom [x_1 x_2 x_3] on a specified element group.
/**
Required input parameters:      

\code{.json}
 "Boundary4":{                  
     "Name": "ExternalForce",   
     "ElementGroup":"top",      
     "dofs": ["pw"],            
     "Force": [1.0e-3]          
},
\endcode

where: \n
'name': The name of this model \n
'ElementGroup': The group of elements this model is being applied to (allows for both boundaries and internal) \n
'dofs': An array of the degrees of freedom that forces are applied to \n
'Force': An array (of the same length as dofs) that indicates the forcing of each degree of freedom \n

*/
class ExternalForceModel: public BaseModel {
    public:
        std::string Name;

        ExternalForceModel(Physics& My_Physics, std::string MyName);
        ~ExternalForceModel();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);
    protected:

    private:
        size_t ElementGroupIndex;
        std::vector<std::string> dofNames;
        std::vector<size_t> dofSteps;
        std::vector<size_t> dofTypes;
        std::vector<double> ForceVals;

};

void Register_ExternalForceModel();
BaseModel* New_ExternalForceModel(Physics& My_Physics, std::string MyNameIn);

#endif

