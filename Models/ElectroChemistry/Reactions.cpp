#include "Reactions.h"

SurfaceReaction::SurfaceReaction(std::string ReactionNameIn, inputData& inputs){
    ReactionName = ReactionNameIn;
    inputs.GetRequired(ReactionType, {"properties","SurfaceReactions",ReactionName,"Type"});
    inputs.GetRequired(SurfaceID, {"properties","SurfaceReactions",ReactionName,"Surface"});

    if (ReactionType == ReactionTypes::ELECTRO_REACTION){
        std::vector<double> Rates; inputs.GetRequired(Rates, {"properties","SurfaceReactions",ReactionName,"i0"});
        i0 = Rates[0];
        i0r = Rates[1];  

        inputs.GetRequired(alpha, {"properties","SurfaceReactions",ReactionName,"alpha"});
        inputs.GetRequired(E_eq, {"properties","SurfaceReactions",ReactionName,"E_eq"});
        inputs.GetRequired(nElectrons, {"properties","SurfaceReactions",ReactionName,"electrons_In"});
    } else {
        throw std::invalid_argument("Undefined reaction type for "+ReactionName+"\n");
    }
    Lumped = false;
    inputs.GetOptional(Lumped, {"properties","SurfaceReactions",ReactionName,"Lumped"});


    std::vector<std::string> InSpeciesNames, OutSpeciesNames;
    std::vector<int> n_in, n_out;
    inputs.GetRequired(AllSpecies, {"properties","Species","Types"});
    inputs.GetRequired(InSpeciesNames, {"properties","SurfaceReactions",ReactionName,"Species_In"});
    inputs.GetRequired(OutSpeciesNames, {"properties","SurfaceReactions",ReactionName,"Species_Out"});
    inputs.GetRequired(n_in, {"properties","SurfaceReactions",ReactionName,"n_In"});
    inputs.GetRequired(n_out, {"properties","SurfaceReactions",ReactionName,"n_Out"});
    inputs.GetRequired(C_ref,  {"properties","SurfaceReactions",ReactionName,"C_ref"});

    nSpecies = AllSpecies.size();
    InSpecies.resize(nSpecies);
    OutSpecies.resize(nSpecies);

    for (size_t i = 0; i < nSpecies; i++){
        InSpecies(i) = 0;
        for (size_t j = 0; j < InSpeciesNames.size(); j++){
            if (InSpeciesNames[j]==AllSpecies[i]){
                InSpecies(i) = n_in[j];
            }
        }
        
        OutSpecies(i) = 0;
        for (size_t j = 0; j < OutSpeciesNames.size(); j++){
            if (OutSpeciesNames[j]==AllSpecies[i]){
                OutSpecies(i) = n_out[j];
            }
        }
    }
}

SurfaceReaction::~SurfaceReaction(){

}

/// @brief 
/// @param SurfaceName input: Name of surface on which reactions occur
/// @param CVec input: Concentration vector
/// @param ePot input: electrolyte potential
/// @param Em input: metal potential
/// @param SpeciesFlux output: reaction flux
/// @param dSpeciesFlux_dC output d(reaction flux)/dC
/// @param dSpeciesFlux_dE output d(reaction flux)/dEm, add a - sign for ePot
/// @param iReact output: reaction current density
/// @param di_dE output: di_dEM, add a - sign for ePot
/// @param di_dC output: di_dC
/// @return does this surface have any reaction
bool SurfaceReaction::GetRates(std::string &SurfaceName, std::vector<double> &CVec, double ePot, double Em, 
                                 std::vector<double> &SpeciesFlux, std::vector<std::vector<double>> &dSpeciesFlux_dC, std::vector<double> &dSpeciesFlux_dE, 
                                 double &iReact, double &di_dE, std::vector<double> &di_dC){
    bool activeSurface = false;

    double rate = 0;
    if (SurfaceName == SurfaceID){
        activeSurface = true;
        double Rate_forward=0.0, Rate_backwards=0.0;
        double kForward=0.0, kBackward=0.0, dk_forward_dePot=0.0, dk_backwards_dePot=0.0;
        double eta=0.0;
        std::vector<double> CVecCopy(nSpecies), dC(nSpecies);

        if (ReactionType == ReactionTypes::ELECTRO_REACTION){
            eta = Em-ePot-E_eq;
            eta = std::max(-1.5, std::min(1.5, eta));
            
            kForward = i0/(F*nElectrons)    * std::exp(-alpha*eta*F/R/T);
            kBackward = i0r/(F*nElectrons) * std::exp((1-alpha)*eta*F/R/T);
            dk_forward_dePot = i0/(F*nElectrons)  * std::exp(-alpha*eta*F/R/T) * -alpha*F/R/T;
            dk_backwards_dePot = i0r/(F*nElectrons) * std::exp((1-alpha)*eta*F/R/T) * (1-alpha)*F/R/T;

            for (size_t i = 0; i < nSpecies; i++){
                CVecCopy[i] = CVec[i];
                dC[i] = 1.0;
                if (CVecCopy[i]<0.0){
                    CVecCopy[i] = 0.0;
                    dC[i] = 0.0;
                }
            }
        }
        Rate_forward = 1.0;
        Rate_backwards = 1.0;
        for (size_t i = 0; i < nSpecies; i++){
            if (InSpecies(i)>0){
                Rate_forward *= CVecCopy[i]/C_ref;
            }
            if (OutSpecies(i)>0){
                Rate_backwards *= CVecCopy[i]/C_ref;
            }
        }
        rate = Rate_forward*kForward-Rate_backwards*kBackward;
        for (size_t i = 0; i < nSpecies; i++){
            SpeciesFlux[i] = rate*(OutSpecies(i)-InSpecies(i));
        }   
        iReact = F*nElectrons*rate;

        // get derivatives
        std::vector<double> dRateForward_dC(nSpecies), dRateBackward_dC(nSpecies);
        
        for (size_t i = 0; i < nSpecies; i++){
            if (InSpecies(i)>0){
                dRateForward_dC[i] = 1.0;
                for (size_t j = 0; j < nSpecies; j++){
                    if (InSpecies(j)>0 && i!=j){
                        dRateForward_dC[i] *= CVecCopy[j]/C_ref;
                    }
                    if (InSpecies(j)>0 && i==j){
                        dRateForward_dC[i] *= 1.0/C_ref;
                    }
                } 
            } else {
                dRateForward_dC[i] = 0.0;
            }
            if (OutSpecies(i)>0){
                dRateBackward_dC[i] = 1.0;
                for (size_t j = 0; j < nSpecies; j++){
                    if (OutSpecies(j)>0 && i!=j){
                        dRateBackward_dC[i] *= CVecCopy[j]/C_ref;
                    }
                    if (OutSpecies(j)>0 && i==j){
                        dRateBackward_dC[i] *= 1.0/C_ref;
                    }
                } 
            } else{
                dRateBackward_dC[i] = 0.0;
            }
        }

        for (size_t i = 0; i < nSpecies; i++){
            dSpeciesFlux_dE[i] = (Rate_forward*dk_forward_dePot-Rate_backwards*dk_backwards_dePot)*(OutSpecies(i)-InSpecies(i));
            di_dC[i] = F*nElectrons*(dRateForward_dC[i]*kForward-dRateBackward_dC[i]*kBackward)*dC[i];
            for (size_t j = 0; j < nSpecies; j++){
                dSpeciesFlux_dC[i][j] = (dRateForward_dC[j]*kForward-dRateBackward_dC[j]*kBackward)*(OutSpecies(i)-InSpecies(i))*dC[j];
            }
        }
        di_dE = F*nElectrons*(Rate_forward*dk_forward_dePot-Rate_backwards*dk_backwards_dePot);
    } else {
        iReact=0.0;
        di_dE = 0.0;
        for (size_t i = 0; i < nSpecies; i++){
            SpeciesFlux[i] = 0.0;
            dSpeciesFlux_dE[i] = 0.0;
            di_dC[i] = 0.0;
            for (size_t j = 0; j < nSpecies; j++){
                dSpeciesFlux_dC[i][j] = 0.0;
            }
        }
    }

    return activeSurface;
}

VolumeReaction::VolumeReaction(std::string ReactionName, inputData& inputs){
    inputs.GetRequired(ReactionType, {"properties","VolumeReactions",ReactionName,"Type"});
    
    if (ReactionType == ReactionTypes::EQUILIBRIUM_REACTION){
        inputs.GetRequired(K, {"properties","VolumeReactions",ReactionName,"K"});
        inputs.GetRequired(k_dummy, {"properties","VolumeReactions",ReactionName,"k_dummy"});
    } else if (ReactionType == ReactionTypes::DYNAMIC_REACTION){
        std::vector<double> Rates; inputs.GetRequired(Rates, {"properties","VolumeReactions",ReactionName,"k"});
        k1 = Rates[0];
        k2 = Rates[1];  
    } else {
        throw std::invalid_argument("Undefined reaction type for "+ReactionName+"\n");
    }

    Lumped = false;
    inputs.GetOptional(Lumped, {"properties","VolumeReactions",ReactionName,"Lumped"});


    std::vector<std::string> InSpeciesNames, OutSpeciesNames;
    std::vector<int> n_in, n_out;
    inputs.GetRequired(AllSpecies, {"properties","Species","Types"});
    inputs.GetRequired(InSpeciesNames, {"properties","VolumeReactions",ReactionName,"Species_In"});
    inputs.GetRequired(OutSpeciesNames, {"properties","VolumeReactions",ReactionName,"Species_Out"});
    inputs.GetRequired(n_in, {"properties","VolumeReactions",ReactionName,"n_In"});
    inputs.GetRequired(n_out, {"properties","VolumeReactions",ReactionName,"n_Out"});
    inputs.GetRequired(C_ref,  {"properties","VolumeReactions",ReactionName,"C_ref"});

    nSpecies = AllSpecies.size();
    InSpecies.resize(nSpecies);
    OutSpecies.resize(nSpecies);

    for (size_t i = 0; i < nSpecies; i++){
        InSpecies(i) = 0;
        for (size_t j = 0; j < InSpeciesNames.size(); j++){
            if (InSpeciesNames[j]==AllSpecies[i]){
                InSpecies(i) = n_in[j];
            }
        }
        
        OutSpecies(i) = 0;
        for (size_t j = 0; j < OutSpeciesNames.size(); j++){
            if (OutSpeciesNames[j]==AllSpecies[i]){
                OutSpecies(i) = n_out[j];
            }
        }
    }
}

VolumeReaction::~VolumeReaction(){

}

bool VolumeReaction::isLumped(){
    return Lumped;
}

double VolumeReaction::GetRates(std::vector<double> &CVec, std::vector<double> &SpeciesFlux, std::vector<std::vector<double>> &dSpeciesFlux, std::vector<size_t> &RelevantSpecies){
    std::vector<double> CVecCopy(CVec.size());
    double rate, Rate_forward=1.0, rate_backwards=1.0;
    double k_forward=0.0, k_backwards=0.0;
    std::vector<double> dC(nSpecies); 
    for (size_t i = 0; i < nSpecies; i++)
    {
        dC[i] = 1.0;
    }

    if (ReactionType == ReactionTypes::EQUILIBRIUM_REACTION){
        k_forward = K*k_dummy;
        k_backwards = k_dummy;

        for (size_t i = 0; i < nSpecies; i++){
            CVecCopy[i] = CVec[i];
            if (CVecCopy[i]<1.0e-16){
                CVecCopy[i] = 1.0e-16;
                dC[i] = 0.0;
            }
        }
    } else if (ReactionType == ReactionTypes::DYNAMIC_REACTION){
        k_forward = k1;
        k_backwards = k2;

        for (size_t i = 0; i < nSpecies; i++){
            CVecCopy[i] = CVec[i];
            if (CVecCopy[i]<0.0){
                CVecCopy[i] = 0.0;
                dC[i] = 0.0;
            }
        }
    }

    Rate_forward = 1.0;
    rate_backwards = 1.0;
    RelevantSpecies.resize(0);
    for (size_t i = 0; i < nSpecies; i++){
        if (InSpecies(i)!=0){
            Rate_forward *= CVecCopy[i]/C_ref;
            RelevantSpecies.push_back(i);
        }
        if (OutSpecies(i)!=0){
            rate_backwards *= CVecCopy[i]/C_ref;
            RelevantSpecies.push_back(i);
        }
    }
    sort( RelevantSpecies.begin(), RelevantSpecies.end() );
    RelevantSpecies.erase( unique( RelevantSpecies.begin(), RelevantSpecies.end() ), RelevantSpecies.end() );

    rate = k_forward*Rate_forward-k_backwards*rate_backwards;
    for (size_t i = 0; i < nSpecies; i++){
        SpeciesFlux[i] = rate*(OutSpecies(i)-InSpecies(i));
    }

    std::vector<double> dRateForward_dC(nSpecies), dRateBackward_dC(nSpecies);
    for (size_t i = 0; i < nSpecies; i++){
        if (InSpecies(i)!=0){
            dRateForward_dC[i] = 1.0;
            for (size_t j = 0; j < nSpecies; j++){
                if (InSpecies(j)!=0 && i!=j){
                    dRateForward_dC[i] *= CVecCopy[j]/C_ref;
                }
                if (InSpecies(j)!=0 && i==j){
                    dRateForward_dC[i] *= 1.0/C_ref;
                }
            } 
        } else {
            dRateForward_dC[i] = 0.0;
        }
        if (OutSpecies(i)!=0){
            dRateBackward_dC[i] = 1.0;
            for (size_t j = 0; j < nSpecies; j++){
                if (OutSpecies(j)!=0 && i!=j){
                    dRateBackward_dC[i] *= CVecCopy[j]/C_ref;
                }
                if (OutSpecies(j)!=0 && i==j){
                    dRateBackward_dC[i] *= 1.0/C_ref;
                }
            } 
        } else{
            dRateBackward_dC[i] = 0.0;
        }
    }

    for (size_t i = 0; i < nSpecies; i++){
        for (size_t j = 0; j < nSpecies; j++){
            dSpeciesFlux[i][j] = (k_forward*dRateForward_dC[j]-k_backwards*dRateBackward_dC[j])*(OutSpecies(i)-InSpecies(i))*dC[j];
            //std::cout << dSpeciesFlux[i][j] << "  ";
        }
        //std::cout << "\n";
    }
    //std::cout << "\n";

    double extraRate = 0;
    std::vector<double> dExtraRatedSpecies(nSpecies);
    
    if (ReactionType == ReactionTypes::EQUILIBRIUM_REACTION){ //protection against severe overshoots
        for (size_t i = 0; i < nSpecies; i++){
            dExtraRatedSpecies[i] = 0.0;
            if (InSpecies(i)>0 && CVec[i]<0){
                extraRate += k_dummy/InSpecies(i)*(1e-6-CVec[i]);
                dExtraRatedSpecies[i] = k_dummy/InSpecies(i)*(-1);
            }
            if (OutSpecies(i)>0 && CVec[i]<0){
                extraRate += -k_dummy/OutSpecies(i)*(1e-6-CVec[i]);
                dExtraRatedSpecies[i] = -k_dummy/OutSpecies(i)*(-1);
            }
        }

        for (size_t i = 0; i < nSpecies; i++){
            SpeciesFlux[i] += extraRate*(OutSpecies(i)-InSpecies(i));
            for (size_t j = 0; j < nSpecies; j++){
                dSpeciesFlux[i][j] += dExtraRatedSpecies[j]*(OutSpecies(i)-InSpecies(i));
            }
        }
    }
    
    return rate;
}