#ifndef ECHEMCONSTRAINT_H
#define ECHEMCONSTRAINT_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"

/// @brief Model that imposes concentration boundary conditions for electro-chemical simulations.
/** 
Required input parameters, within the "Models" node:

\code{.json}
    "Boundary1":{
        "Name": "Electrochemistry/EchemConstraints",
        "NodeGroup_C":"Left",
        "NodeGroup_E":"Left",
        "Species": ["H", "OH", "Fe", "FeOH", "Na", "Cl", "O2"],
        "ePot": 0.0
    },
\endcode

and within the material parameters node, for each species listed:

\code{.json}
    "Cl":{
        "D": 2.0e-9,
        "z": -1,
        "C0": 10.0,
        "Aqueous": true
    },
    "O2":{
        "D": 1.0e-9,
        "z": 0,
        "C0": 1.0,
        "Aqueous": false
    }
\endcode

Where the concentrations "C0" are used as imposed boundary concentration. One species can be defined as being resolved through electroneutrality using ""C0": "Neutrality"". This DOES NOT alter how this species is simulated, it merely calculates the imposed concentration based on the boundary conditions detailed in the file.
*/
class EchemConstraints: public BaseModel{
    public:
        std::vector<std::string> DofNames_C;
        std::string DofNames_ePot = "ePot";

        EchemConstraints(Physics& My_Physics, std::string MyName);
        ~EchemConstraints();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);
        void Commit(int CommitType);
    protected:

    private:
        size_t nSpecies;
        std::vector<std::string> Species;
        bool ChemPotBased;

        std::vector<double> z, C0;
        double E0;

        size_t NodeGroupIndex_C, NodeGroupIndex_E;
        std::vector<size_t> dofTypes_C;
        size_t dofType_ePot;
        size_t Step_C;

        double F = 96485.3321;
        double R = 8.31446261815324;
        double T = 293.15;
};

void Register_EchemConstraints();
BaseModel* New_EchemConstraints(Physics& My_Physics, std::string MyNameIn);

#endif

