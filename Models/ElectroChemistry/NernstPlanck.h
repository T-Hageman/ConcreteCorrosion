#ifndef NERNSTPLANCK_H
#define NERNSTPLANCK_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"
#include "Reactions.h"

/// @brief Model that implements the Nernst-Planck equations for multiple species.
/** 
Required input parameters, within the "Models" node:

\code{.json}
    "NernstPlanck":{
        "Name": "Electrochemistry/NernstPlanck",
        "ElementGroup_C": "Interior",
        "ElementGroup_E": "Interior",
        "Porosity": 0.001,
        "PorosityFactor": 1.5,
        "Sw": 0.25,
        "S_irr": 0.2,
        "SatFactor": 2.0,
        "Stabilisation": "None"
    },
\endcode
where "ElementGroup_C" and "ElementGroup_E" refer to the name of the elements that define the domain, "Porosity" the porosity of the concrete (set to 1 to simulate a free-flowing electrolyte), "PorosityFactor" the exponent used to determine the effective diffusivity, 
"Sw" and "S_irr" the current saturation and saturation below which pores become disconnected, and "Stabilisation" indicating which stabilisation scheme will be used.


Within the properties node, for each species listed, define the diffusivity, particle charge, initial concentration, and whether the species is aqueous (assumes gasseous otherwise):

\code{.json}
    "Species":{
        "Types": ["H", "OH", "Fe", "FeOH", "Na", "Cl", "O2"],
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
    }
\endcode

The relations solved by this model are the Nernst-Planck equation:
\f[ S^*_\pi \phi \dot{C}_\pi +\mathbf{\nabla}\cdot\left(-D_\pi^{\text{eff}} \mathbf{\nabla}C_\pi\right)+\frac{z_\pi F}{RT}\mathbf{\nabla}\cdot\left(-D_\pi^{\text{eff}} C_\pi \mathbf{\nabla}\varphi\right) + \phi S^* R_\pi = 0   \f]
with the effective diffusivity given by:
\f[ D_\pi^{\text{eff}} = \frac{D_\pi \phi}{\tau} S_\text{eff}^2 = \phi^{3/2} D_\pi \left(\frac{S_\text{w}-S_\text{irr}}{1-S_\text{irr}}\right)^2  \f]

Additionally, this model considers the electroneutrality condition as:
\f[ \sum_\pi \phi z_\pi C_\pi = 0  \f]

Finally, this model implements volume reactions taking place within the pore, see the description of VolumeReaction below for more details.
*/
class NernstPlanck: public BaseModel{
    public:
        std::vector<std::string> DofNames_C;
        std::string DofNames_ePot = "ePot";

        NernstPlanck(Physics& My_Physics, std::string MyName);
        ~NernstPlanck();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);
        void Commit(int CommitType);
    protected:

    private:
        size_t nSpecies, nReactions;
        bool ChemPotBased;
        std::vector<std::string> Species;
        std::vector<VolumeReaction> Reactions;

        double porosity, PorosityFactor;
        double Sw, S_irr, S_factor;
        std::vector<double> z;
        std::vector<double> capacity;
        std::vector<Eigen::Matrix3d> D;

        size_t ElemGroupIndex_C, ElemGroupIndex_E;
        std::vector<size_t> dofTypes_C;
        size_t dofType_ePot;
        size_t Step_C;

        double F = 96485.3321;
        double R = 8.31446261815324;
        double T = 293.15;
};

void Register_NernstPlanck();
BaseModel* New_NernstPlanck(Physics& My_Physics, std::string MyNameIn);

#endif

