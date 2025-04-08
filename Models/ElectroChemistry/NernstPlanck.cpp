#include "NernstPlanck.h"
#include "../../Physics/physics.h"

void Register_NernstPlanck(){
    ModelNames.push_back("Electrochemistry/NernstPlanck");
    ModelCreators.push_back(New_NernstPlanck);
}
BaseModel* New_NernstPlanck(Physics& My_Physics, std::string MyNameIn){
    return new NernstPlanck(My_Physics, MyNameIn);
}

NernstPlanck::NernstPlanck(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "Electrochemistry/NernstPlanck";
}

NernstPlanck::~NernstPlanck(){

}

void NernstPlanck::init(inputData& inputs){
    Setup(inputs);

    std::vector<size_t> UniqueNodes_C = mesh->ElementGroups[ElemGroupIndex_C].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes_C, dofTypes_C);

    std::vector<size_t> UniqueNodes_E = mesh->ElementGroups[ElemGroupIndex_E].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes_E, dofType_ePot);
}

void NernstPlanck::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
}

void NernstPlanck::save(SaveDataFile& data){

}

void NernstPlanck::Commit(int CommitType){

}

void NernstPlanck::Setup(inputData& inputs){
    std::string groupName_C, groupName_E;
    inputs.GetRequired(groupName_C, {"Models",MyName,"ElementGroup_C"});
    inputs.GetRequired(groupName_E, {"Models",MyName,"ElementGroup_E"});

    ElemGroupIndex_C = mesh->GetElementGroupIdx(groupName_C);
    ElemGroupIndex_E = mesh->GetElementGroupIdx(groupName_E);

    bool hasPorosity = inputs.GetOptional(porosity, {"Models",MyName,"Porosity"});
    if (hasPorosity){
        inputs.GetRequired(PorosityFactor, {"Models",MyName,"PorosityFactor"});
    } else {
        porosity = 1.0;
        PorosityFactor = 1.0;
    }

    bool hasSaturation = inputs.GetOptional(Sw, {"Models",MyName,"Sw"});
    if (hasSaturation){
        inputs.GetRequired(S_irr, {"Models",MyName,"S_irr"});
        inputs.GetRequired(S_factor, {"Models",MyName,"SatFactor"});
        if (Sw<S_irr){
            Sw = S_irr;
        }
    } else {
        Sw = 1.0;
        S_irr = 0.0;
        S_factor = 1.0;
    }

    ChemPotBased = false;
    inputs.GetOptional(ChemPotBased, {"properties","Species","UseChemPotentials"});

    inputs.GetRequired(Species, {"properties","Species","Types"});
    nSpecies = Species.size();

    D.resize(nSpecies);
    z.resize(nSpecies);
    capacity.resize(nSpecies);
    DofNames_C.resize(nSpecies);
    double D_scalar;
    bool aqueous;
    for (size_t i = 0; i < nSpecies; i++){
        DofNames_C[i] = Species[i];
        inputs.GetRequired(D_scalar,{"properties","Species",Species[i],"D"});
        aqueous = true;
        inputs.GetOptional(aqueous,{"properties","Species",Species[i],"Aqueous"});
        capacity[i] = 1.0;

        D_scalar *= std::pow(porosity,PorosityFactor);
        capacity[i] *= porosity;
        if (aqueous){
            D_scalar *= std::pow((Sw-S_irr)/(1.0-S_irr), S_factor);
            capacity[i] *= Sw;
        }
        D[i] = D_scalar * Eigen::Matrix3d::Identity();
        inputs.GetRequired(z[i], {"properties","Species",Species[i],"z"});
    }

    std::vector<std::string> ReactionNames; inputs.GetRequired(ReactionNames, {"properties","VolumeReactions","Reactions"});

    nReactions = ReactionNames.size();
    for (size_t i = 0; i < nReactions; i++){
        Reactions.push_back(VolumeReaction(ReactionNames[i], inputs));
    }


    std::vector<size_t> dofSteps;
    dofs->getDofTypesSteps(DofNames_C, dofTypes_C, dofSteps);
    Step_C = dofSteps[0];
    for (size_t i = 1; i < nSpecies; i++){
        if (Step_C!=dofSteps[i]){
            throw std::invalid_argument(ModelName+" requires "+ DofNames_C[0] + " and " + DofNames_C[i] + " to be in the same solver step\n");
        }
    }
    size_t EStep;
    dofs->getDofTypesSteps(DofNames_ePot, dofType_ePot, EStep);
    if (EStep != Step_C){
        throw std::invalid_argument(ModelName+" requires "+ DofNames_C[0] + " and " + DofNames_ePot + " to be in the same solver step\n");
    }
}

void NernstPlanck::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    //double t = physics->time;
    if (step == Step_C){
        size_t ipcount = mesh->ElementGroups[ElemGroupIndex_C].BaseElem->ipcount;
        size_t nNodes_C = mesh->ElementGroups[ElemGroupIndex_C].NNodes_per_elem;
        size_t nNodes_E = mesh->ElementGroups[ElemGroupIndex_E].NNodes_per_elem;

        std::vector<double> w(ipcount);

        std::vector<Eigen::RowVectorXd> Nc(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nc[ip].resize(nNodes_C); 
        std::vector<Eigen::MatrixXd> Gc(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gc[ip].resize(mesh->dim, nNodes_C); 

        std::vector<Eigen::RowVectorXd> Ne(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Ne[ip].resize(nNodes_E); 
        std::vector<Eigen::MatrixXd> Ge(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Ge[ip].resize(mesh->dim, nNodes_E); 

        Eigen::VectorXd WLumped_C(nNodes_C), WLumped_E(nNodes_E);

        std::vector<size_t> El_C(nNodes_C), El_E(nNodes_E);

        std::vector<PetscInt> dofsE(nNodes_E);
        std::vector<std::vector<PetscInt>> dofsC(nSpecies); for (size_t s=0; s<nSpecies; s++) dofsC[s].resize(nNodes_C);

        Eigen::VectorXd E(nNodes_E), E_Old(nNodes_E);
        std::vector<Eigen::VectorXd> C(nSpecies), COld(nSpecies);
        for (size_t i = 0; i < nSpecies; i++){
            C[i].resize(nNodes_C);
            COld[i].resize(nNodes_C);
        }

        std::vector<double> cLoc(nSpecies), dcLoc_dt(nSpecies);
        std::vector<Eigen::VectorXd> gradC(nSpecies); for (size_t s=0; s<nSpecies; s++) gradC[s].resize(mesh->dim);
        std::vector<Eigen::RowVectorXd> dC_dDof(nSpecies); for (size_t s=0; s<nSpecies; s++) dC_dDof[s].resize(nNodes_C);
        std::vector<Eigen::MatrixXd> dC_dGradDof(nSpecies); for (size_t s=0; s<nSpecies; s++) dC_dGradDof[s].resize(mesh->dim, nNodes_C);


        std::vector<double> ReactionFlux(nSpecies);
        std::vector<std::vector<double>> dReactionFlux(nSpecies); for (size_t i = 0; i < nSpecies; i++) dReactionFlux[i].resize(nSpecies);
        std::vector<size_t> RelevantSpecies;
        
        std::vector<Eigen::VectorXd> F_C(nSpecies); for (size_t s=0; s<nSpecies; s++) F_C[s].resize(nNodes_C);
        Eigen::VectorXd F_E(nNodes_E);
        Eigen::MatrixXd K_EE(nNodes_E, nNodes_E);
        std::vector<Eigen::MatrixXd> K_EC(nSpecies); for (size_t s=0; s<nSpecies; s++) K_EC[s].resize(nNodes_E, nNodes_C);
        std::vector<Eigen::MatrixXd> K_CE(nSpecies); for (size_t s=0; s<nSpecies; s++) K_CE[s].resize(nNodes_C, nNodes_E);
        std::vector<std::vector<Eigen::MatrixXd>> K_CC(nSpecies);
        for (size_t i = 0; i < nSpecies; i++){
            K_CC[i].resize(nSpecies);
            for (size_t j = 0; j < nSpecies; j++){
                K_CC[i][j].resize(nNodes_C, nNodes_C);
            }
        }

        double dt = physics->timeScheme->dt;

        for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndex_C].NElems; el++){
            mesh->getShapeGrads(ElemGroupIndex_C, el, El_C, Nc, Gc, w);
            mesh->getShapeGrads(ElemGroupIndex_E, el, El_E, Ne, Ge, w);

            dofs->getDofForNodes(El_E, dofType_ePot, dofsE);
            physics->StateVectors[Step_C].GetValues(dofsE, E);
            physics->StateVectorsOld[Step_C].GetValues(dofsE, E_Old);

            for (size_t s = 0; s < nSpecies; s++){
                dofs->getDofForNodes(El_C, dofTypes_C[s], dofsC[s]);
                physics->StateVectors[Step_C].GetValues(dofsC[s], C[s]);    
                physics->StateVectorsOld[Step_C].GetValues(dofsC[s], COld[s]);  
            }

            F_E.setZero();
            K_EE.setZero();
            for (size_t s = 0; s < nSpecies; s++){
                F_C[s].setZero();
                K_EC[s].setZero();
                K_CE[s].setZero();
                for (size_t s2 = 0; s2 < nSpecies; s2++){
                    K_CC[s][s2].setZero();
                }
            }
            WLumped_C.setZero();
            WLumped_E.setZero();

            for (size_t ip = 0; ip < ipcount; ip++){

                for (size_t s = 0; s < nSpecies; s++){
                    if (ChemPotBased){
                        cLoc[s] = std::exp(Nc[ip]*C[s]);
                        dcLoc_dt[s] = std::exp(Nc[ip]*C[s])-std::exp(Nc[ip]*COld[s]);
                        dC_dDof[s] = std::max(1.0e-4,std::exp(Nc[ip]*C[s]))*Nc[ip];
                        gradC[s] = std::exp(Nc[ip]*C[s])*Gc[ip]*C[s];
                        dC_dGradDof[s] = std::exp(Nc[ip]*C[s])*Gc[ip] + std::exp(Nc[ip]*C[s])*Gc[ip]*C[s]*Nc[ip];
                    } else {
                        cLoc[s] = Nc[ip]*C[s];
                        dcLoc_dt[s] = Nc[ip]*(C[s]-COld[s]);
                        dC_dDof[s] = Nc[ip];
                        gradC[s] = Gc[ip]*C[s];
                        dC_dGradDof[s] = Gc[ip];
                    }
                    
                    Eigen::Matrix3d Ds = D[s];

                    //Diffusion
                    F_C[s]      += w[ip]*Gc[ip].transpose()*Ds*gradC[s];
                    K_CC[s][s]  += w[ip]*Gc[ip].transpose()*Ds*dC_dGradDof[s];

                    //electromigration
                    double zF_RT = z[s]*F/R/T;
                    F_C[s]      += w[ip]*Gc[ip].transpose()*D[s]*zF_RT*(cLoc[s])*Ge[ip]*E;
                    K_CC[s][s]  += w[ip]*Gc[ip].transpose()*D[s]*zF_RT*(Ge[ip]*E)*dC_dDof[s];
                    K_CE[s]     += w[ip]*Gc[ip].transpose()*D[s]*zF_RT*(cLoc[s])*Ge[ip];
                }
                
                //reactions
                for (size_t r = 0; r < Reactions.size(); r++){
                    if (Reactions[r].isLumped() == false){
                        Reactions[r].GetRates(cLoc, ReactionFlux, dReactionFlux, RelevantSpecies);
                        for (size_t s = 0; s < RelevantSpecies.size(); s++){
                            size_t sa = RelevantSpecies[s];
                            F_C[sa] -= w[ip]*Nc[ip].transpose()*ReactionFlux[sa]*capacity[sa];
                            for (size_t s2 = 0; s2 < nSpecies; s2++){
                                size_t sb = RelevantSpecies[s2];
                                K_CC[sa][sb] -= w[ip]*Nc[ip].transpose()*dC_dDof[sb]*dReactionFlux[sa][sb]*capacity[sa];
                            }
                        }
                    }
                }
                
                WLumped_C += w[ip]*Nc[ip];
                WLumped_E += w[ip]*Ne[ip];
            }

            //lumped integrations
            for (size_t n = 0; n < nNodes_C; n++){
                for (size_t s = 0; s < nSpecies; s++){
                    if (ChemPotBased){
                        cLoc[s] = std::exp(C[s](n));
                        dcLoc_dt[s] = std::exp(C[s](n))-std::exp(COld[s](n));
                        dC_dDof[s](n) = std::max(1.0e-4,std::exp(C[s](n)));
                    } else {
                        cLoc[s] = C[s](n);
                        dcLoc_dt[s] = C[s](n)-COld[s](n);
                        dC_dDof[s](n) = 1.0;
                    }
                }

                //Capacity
                for (size_t s = 0; s < nSpecies; s++){
                    F_C[s](n)       += WLumped_C[n]*capacity[s]*dcLoc_dt[s]/dt;
                    K_CC[s][s](n,n) += WLumped_C[n]*capacity[s]*dC_dDof[s](n)/dt;
                }

                //reactions
                for (size_t r = 0; r < Reactions.size(); r++){
                    if (Reactions[r].isLumped() == true){
                        Reactions[r].GetRates(cLoc, ReactionFlux, dReactionFlux, RelevantSpecies);
                        for (size_t s = 0; s < RelevantSpecies.size(); s++){
                            size_t sa = RelevantSpecies[s];
                            F_C[sa](n) -= WLumped_C[n]*porosity*ReactionFlux[sa];
                            for (size_t s2 = 0; s2 < RelevantSpecies.size(); s2++){
                                size_t sb = RelevantSpecies[s2];
                                K_CC[sa][sb](n,n) -= WLumped_C[n]*porosity*dC_dDof[sb](n)*dReactionFlux[sa][sb];
                            }
                        }
                    }
                }
            }

            //electroneutrality
            for (size_t n = 0; n < nNodes_E; n++){
                for (size_t s = 0; s < nSpecies; s++){
                    if (C[s](n)>=0){
                        F_E(n)         += WLumped_E[n]*z[s]*C[s](n);
                        K_EC[s](n,n)   += WLumped_E[n]*z[s];
                    }
                }
            }            

            //adding to matrix
            VecAdd(f, dofsE, F_E);
            MatAdd(K, dofsE, dofsE, K_EE);

            for (size_t i = 0; i < nSpecies; i++){
                MatAdd(K, dofsE, dofsC[i], K_EC[i]);

                VecAdd(f, dofsC[i], F_C[i]);
                MatAdd(K, dofsC[i], dofsE, K_CE[i]);
                for (size_t j = 0; j < nSpecies; j++){
                    if (K_CC[i][j].isZero(1.0e-30) == false){
                        MatAdd(K, dofsC[i], dofsC[j], K_CC[i][j]);
                    }
                }
            }
        }
    }
}