#ifndef ELECTROSURFACE_H
#define ELECTROSURFACE_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"
#include "Reactions.h"

/// @brief This model adds the reactions taking place on the metal surface, and includes the current-conservation condition.
/** 
Required input parameters, within the "Models" node:

\code{.json}
        "Surface":{
            "Name": "Electrochemistry/ElectroSurface",
            "AreaGroups_C": ["Pit", "Bar"],
            "AreaGroups_E": ["Pit", "Bar"],
            "SurfaceReactions": [["Cathode","Anode"],["Cathode"]],
            "ChargeConservation": true,
            "E_m": 0.0,
            "ActiveCurrent_Threshold": 1.0e-5
        },
\endcode
where "AreaGroups_C" and "AreaGroups_E" provide the names of the elements used to define the domain (can be the same for both), and SurfaceReactions details which reactions are allowed to occur on which group.
The boolean "ChargeConservation" indicates whther a metal potential is directly imposed, or if the current-conservation condition is imposed as:
\f[ \sum I_\pi+I_\text{ext}=\int_\Gamma \sum i_\pi\;\text{d}\Gamma +I_\text{ext}  = 0 \f]
here, assuming no external current.

For the definition of the surface reactions, see the description of SurfaceReaction below for more details.
*/
class ElectroSurface: public BaseModel{
    public:
        std::vector<std::string> DofNames_C;
        std::string DofNames_ePot = "ePot";
        std::string DofName_EM = "Em";

        ElectroSurface(Physics& My_Physics, std::string MyName);
        ~ElectroSurface();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);
        void Commit(int CommitType);

        size_t hasTimeData(std::vector<std::string>& DataNames);
        void GetTimeData(std::vector<double>& DataValues);

        bool SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data);
    protected:

    private:
        size_t nSpecies, nReactions;
        std::vector<std::string> Species;

        std::vector<size_t> ElemGroupIndices_C, ElemGroupIndices_E;
        std::vector<std::vector<std::string>> InterfaceReactionTypes;
        std::vector<size_t> dofTypes_C;
        size_t dofType_ePot;
        size_t Step_C;

        size_t nAreas; 
        std::vector<SurfaceReaction> Reactions;
        std::vector<double> ReactionCurrents;
        std::vector<double> ReactionArea;
        double ActiveCurrent_Threshold;


        size_t EMNode = 0;
        size_t dofType_EM;
        double E_m;
        bool ChargeConservation;
        bool ChemPotBased;

        double F = 96485.3321;
        double R = 8.31446261815324;
        double T = 293.15;
        double Curr_Scale = 1.0e8;
};

void Register_ElectroSurface();
BaseModel* New_ElectroSurface(Physics& My_Physics, std::string MyNameIn);

#endif

