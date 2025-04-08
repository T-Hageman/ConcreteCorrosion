#ifndef REACTIONS_H
#define REACTIONS_H

#include "../../InputsOutputs/inputData.h"
#include "../../InputsOutputs/SaveData.h"

/// @brief Indicator of type of reaction being considered
namespace ReactionTypes{
    /// @brief Reaction is dynamic, considering forwards and backwards rate
    constexpr std::string_view DYNAMIC_REACTION = "Dynamic";

    /// @brief Reaction is in equilibrium, considering equilibrium constant and dummy rate
    constexpr std::string_view EQUILIBRIUM_REACTION = "Equilibrium";

    /// @brief Electro-chemical reaction, only valid for surfaces
    constexpr std::string_view ELECTRO_REACTION = "Electrochemical";
}

/// @brief Evaluates reaction rates for surface reactions
/** 
Each reaction is defined within the parameter node as

\code{.json}
    "Reactions": ["HER","OER","Corrosion"],
    "HER":{
        "Reaction": "2H+ + 2e- <-> H2",
        "Type": "Electrochemical",
        "i0": [1.0e-2, 0.0],
        "alpha": 0.5,
        "E_eq": 0.0,
        "electrons_In": 2,
        "Species_In":["H"],
        "n_In":[2],
        "Species_Out":[],
        "n_Out":[],
        "C_ref": 1.0e3,
        "Surface": "Cathode"
    }, 
\endcode
where:
- Reaction is not used, and is treated as a comment field
- Type indicates which type of reaction this is (and which parameters it will check for)
- i0 is the reaction rate constant, given as a current for the forward and backward rate, each can be set to zero for uni-directional reactions
- alpha is the charge transfer coefficient
- E_eq is the Equilibrium potential
- Electrons_in indicates how many electrons are involved (positive for being on the input side)
- Species_in lists the species that are used as input for this reaction
- n_In lists the amount of species involved
- Species_out defines the output species
- n_Out defines the amount of outputs for each species
- C_ref is the reference concentration used, typically 1000 mol/m3 (1 M)
- Surface indicates the surface this reaction takes place on, with these surfaces being defined in ElectroSurface

Each reaction rate is evaluated as:
\f[  \nu =  \frac{i_0(0) e^{-\alpha\frac{\eta F}{RT}}}{n_{electrons}F} \prod C_{in}/C_{ref}  -   \frac{i_0(1) e^{(1-\alpha)\frac{\eta F}{RT}} }{n_{electrons}F} \prod C_{out}/C_{ref}  \f]
with the resultant species fluxes then given by:
\f[ \nu_\pi = (n_{out} - n_{in} ) \nu  \f]
*/
class SurfaceReaction{
    public:
        SurfaceReaction(std::string ReactionName, inputData& inputs);
        ~SurfaceReaction();
        bool GetRates(std::string &SurfaceName, std::vector<double> &CVec, double ePot, double Em, 
                        std::vector<double> &SpeciesFlux, std::vector<std::vector<double>> &dSpeciesFlux_dC, std::vector<double> &dSpeciesFlux_dE, 
                        double &i, double &di_dE, std::vector<double> &di_dC);

        std::string ReactionName;
    private:
        std::vector<std::string> AllSpecies;
        Eigen::VectorXi InSpecies, OutSpecies;
        size_t nSpecies;

        std::string ReactionType;
        bool Lumped;

        double i0, i0r;
        double C_ref;
        double alpha, E_eq;
        int nElectrons;
        std::string SurfaceID;  

        double F = 96485.3321;
        double R = 8.31446261815324;
        double T = 293.15;  
};


/// @brief Evaluates reaction rates for volume reactions
/** 
Each reaction is defined within the parameter node as

\code{.json}
    "VolumeReactions":{
        "Reactions": ["Auto-Ionisation","Corrosion1","Corrosion2"],
        "Auto-Ionisation":{
            "Type": "Equilibrium",
            "K": 1.0e-14,
            "k_dummy": 1.0e7,
            "Species_In":[],
            "n_In":[],
            "Species_Out":["H","OH"],
            "n_Out":[1,1],
            "C_ref": 1.0e3,
            "Lumped": true
        },
        "Corrosion1":{
            "Type": "Dynamic",
            "k": [1.0e1, 1.0e1],
            "Species_In": ["Fe"],
            "n_In": [1],
            "Species_Out": ["FeOH","H"],
            "n_Out": [1,1],
            "C_ref": 1.0e3,
            "Lumped": true
        }
    }
\endcode
where:
- Reactions lists all reactions taking place
- Type indicates which type of reaction this is (and which parameters it will check for), Equilibrium or Dynamic
- Species_in lists the species that are used as input for this reaction
- n_In lists the amount of species involved
- Species_out defines the output species
- n_Out defines the amount of outputs for each species
- C_ref is the reference concentration used, typically 1000 mol/m3 (1 M)
- Lumped indicates whether lumped integration should be used to stabilise fast reactions

In case of Equilibrium reactions, aditionally required:
- K, the equilbrium constant
- k_dummy, the dummy parameter, chosen high enough to enforce equilibrium

For dynamic reactions, define:
- k, a vector indicating the forward and reverse reaction rates
*/
class VolumeReaction{
    public:
        VolumeReaction(std::string ReactionName, inputData& inputs);
        ~VolumeReaction();
        double GetRates(std::vector<double> &CVec, std::vector<double> &SpeciesFlux, std::vector<std::vector<double>> &dSpeciesFlux, std::vector<size_t> &RelevantSpecies);
        bool isLumped();

    private:
        std::vector<std::string> AllSpecies;
        Eigen::VectorXi InSpecies, OutSpecies;
        size_t nSpecies;

        std::string ReactionType;
        bool Lumped;
        double K, k_dummy;
        double k1, k2;
        double C_ref;
};

#endif