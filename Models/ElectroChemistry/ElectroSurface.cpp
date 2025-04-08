#include "ElectroSurface.h"
#include "../../Physics/physics.h"

void Register_ElectroSurface(){
    ModelNames.push_back("Electrochemistry/ElectroSurface");
    ModelCreators.push_back(New_ElectroSurface);
}
BaseModel* New_ElectroSurface(Physics& My_Physics, std::string MyNameIn){
    return new ElectroSurface(My_Physics, MyNameIn);
}

ElectroSurface::ElectroSurface(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "Electrochemistry/ElectroSurface";
}

ElectroSurface::~ElectroSurface(){

}

void ElectroSurface::init(inputData& inputs){
    Setup(inputs);

    for (size_t i = 0; i < nAreas; i++){
        std::vector<size_t> UniqueNodes_C = mesh->ElementGroups[ElemGroupIndices_C[i]].GetUniqueNodes();
        dofs->AddDofs(UniqueNodes_C, dofTypes_C);

        std::vector<size_t> UniqueNodes_E = mesh->ElementGroups[ElemGroupIndices_E[i]].GetUniqueNodes();
        dofs->AddDofs(UniqueNodes_E, dofType_ePot);
    }

    if (ChargeConservation){
        dofs->AddDofs(EMNode, dofType_EM);
    }
}

void ElectroSurface::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
}

void ElectroSurface::save(SaveDataFile& data){

}

void ElectroSurface::Commit(int CommitType){

}

void ElectroSurface::Setup(inputData& inputs){
    std::vector<std::string> SurfNames_C, SurfNames_E;
    inputs.GetRequired(SurfNames_C, {"Models",MyName,"AreaGroups_C"});
    inputs.GetRequired(SurfNames_E, {"Models",MyName,"AreaGroups_E"});
    inputs.GetRequired(InterfaceReactionTypes, {"Models",MyName,"SurfaceReactions"});

    nAreas = SurfNames_C.size();
    ElemGroupIndices_C.resize(nAreas);
    ElemGroupIndices_E.resize(nAreas);
    for (size_t i = 0; i < nAreas; i++){
        ElemGroupIndices_C[i] = mesh->GetElementGroupIdx(SurfNames_C[i]);
        ElemGroupIndices_E[i] = mesh->GetElementGroupIdx(SurfNames_E[i]);
    }

    inputs.GetRequired(Species, {"properties","Species","Types"});
    nSpecies = Species.size();
    DofNames_C.resize(nSpecies);
    for (size_t i = 0; i < nSpecies; i++){
        DofNames_C[i] = Species[i];
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

    ChargeConservation = false;
    inputs.GetOptional(ChargeConservation, {"Models",MyName,"ChargeConservation"});
    if (ChargeConservation){
        dofs->getDofTypesSteps(DofName_EM, dofType_EM, EStep);
        if (EStep != Step_C){
            throw std::invalid_argument(ModelName+" requires "+ DofNames_C[0] + " and " + DofName_EM + " to be in the same solver step\n");
        }
    } else {
        inputs.GetRequired(E_m, {"Models",MyName,"E_m"});
    }

    std::vector<std::string> ReactionNames; inputs.GetRequired(ReactionNames, {"properties","SurfaceReactions","Reactions"});
    nReactions = ReactionNames.size();
    ReactionCurrents.resize(nReactions); 
    ReactionArea.resize(nReactions);
    for (size_t i = 0; i < nReactions; i++){
        Reactions.push_back(SurfaceReaction(ReactionNames[i], inputs));
    }
    ActiveCurrent_Threshold = 1.0e-5;
    inputs.GetOptional(ActiveCurrent_Threshold, {"Models",MyName,"ActiveCurrent_Threshold"});

    ChemPotBased = false;
    inputs.GetOptional(ChemPotBased, {"properties","Species","UseChemPotentials"});
}

void ElectroSurface::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    //double t = physics->time;
    if (step == Step_C){
        size_t Egroup_C, Egroup_E;
        size_t nReactionTypes;
        std::vector<std::string> SurfaceTypes;

        size_t ipcount = mesh->ElementGroups[ElemGroupIndices_C[0]].BaseElem->ipcount;
        size_t nNodes_C = mesh->ElementGroups[ElemGroupIndices_C[0]].NNodes_per_elem;
        size_t nNodes_E = mesh->ElementGroups[ElemGroupIndices_E[0]].NNodes_per_elem;

        std::vector<double> w(ipcount);

        std::vector<Eigen::RowVectorXd> Nc(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Nc[ip].resize(nNodes_C); 
        std::vector<Eigen::MatrixXd> Gc(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Gc[ip].resize(mesh->dim-1, nNodes_C); 

        std::vector<Eigen::RowVectorXd> Ne(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Ne[ip].resize(nNodes_E); 
        std::vector<Eigen::MatrixXd> Ge(ipcount); for (size_t ip=0; ip<ipcount;  ip++) Ge[ip].resize(mesh->dim-1, nNodes_E); 

        Eigen::VectorXd WLumped(nNodes_C);

        std::vector<size_t> El_C(nNodes_C), El_E(nNodes_E);
        std::vector<PetscInt> dofsE(nNodes_E);
        std::vector<std::vector<PetscInt>> dofsC(nSpecies); for (size_t s=0; s<nSpecies; s++) dofsC[s].resize(nNodes_C);
        
        Eigen::VectorXd E(nNodes_E);
        std::vector<Eigen::VectorXd> C(nSpecies);
        for (size_t i = 0; i < nSpecies; i++){
            C[i].resize(nNodes_C);
        }

        std::vector<double> cLoc(nSpecies);
        std::vector<Eigen::RowVectorXd> dC_dDof(nSpecies); for (size_t s=0; s<nSpecies; s++) dC_dDof[s].resize(nNodes_C);

        std::vector<double> ReactionFlux(nSpecies);
        std::vector<std::vector<double>> dReactionFlux_dC(nSpecies); for (size_t i = 0; i < nSpecies; i++) dReactionFlux_dC[i].resize(nSpecies);
        std::vector<double> dReactionFlux_dE(nSpecies);
        double iReact, diReact_dE; 
        std::vector<double> diReact_dC(nSpecies);
        std::vector<size_t> RelevantSpecies;
        double eLoc;
        bool hasReaction;

        std::vector<Eigen::VectorXd> F_C(nSpecies); for (size_t s=0; s<nSpecies; s++) F_C[s].resize(nNodes_C);
        std::vector<std::vector<Eigen::MatrixXd>> K_CC(nSpecies);
        for (size_t i = 0; i < nSpecies; i++){
            K_CC[i].resize(nSpecies);
            for (size_t j = 0; j < nSpecies; j++){
                K_CC[i][j].resize(nNodes_C, nNodes_C);
            }
        }
        std::vector<Eigen::MatrixXd> K_CE(nSpecies); for (size_t s=0; s<nSpecies; s++) K_CE[s].resize(nNodes_C, nNodes_E);
        std::vector<Eigen::MatrixXd> K_CEM(nSpecies); for (size_t s=0; s<nSpecies; s++) K_CEM[s].resize(nNodes_C, 1);

        double F_EM;
        double K_EMEM;
        Eigen::MatrixXd K_EME(1,nNodes_E);
        std::vector<Eigen::MatrixXd> K_EMC(nSpecies); for (size_t s=0; s<nSpecies; s++) K_EMC[s].resize(1, nNodes_C);

        PetscInt dofEM;
        double EM;
        if (ChargeConservation){
            dofs->getDofForNodes(EMNode, dofType_EM, dofEM);
            EM = physics->StateVectors[Step_C].GetValue(dofEM);
        } else {
            EM = E_m;
        }
        
        for (size_t i = 0; i < nReactions; i++){
            ReactionCurrents[i] = 0.0;
            ReactionArea[i] = 0.0;
        }

        for (size_t area = 0; area < nAreas; area++){
            Egroup_C = ElemGroupIndices_C[area];
            Egroup_E = ElemGroupIndices_E[area];
            nReactionTypes = InterfaceReactionTypes[area].size();
            SurfaceTypes.resize(nReactionTypes);
            for (size_t i = 0; i < nReactionTypes; i++){
                SurfaceTypes[i] = InterfaceReactionTypes[area][i];
            }
            
            for (size_t el = 0; el < mesh->ElementGroups[Egroup_C].NElems; el++){
                mesh->getShapeGrads(Egroup_C, el, El_C, Nc, Gc, w);
                mesh->getShapeGrads(Egroup_E, el, El_E, Ne, Ge, w);

                dofs->getDofForNodes(El_E, dofType_ePot, dofsE);
                physics->StateVectors[Step_C].GetValues(dofsE, E);

                for (size_t s = 0; s < nSpecies; s++){
                    dofs->getDofForNodes(El_C, dofTypes_C[s], dofsC[s]);
                    physics->StateVectors[Step_C].GetValues(dofsC[s], C[s]);    
                }

                F_EM = 0.0;
                K_EMEM = 0.0;
                K_EME.setZero();
                for (size_t s = 0; s < nSpecies; s++){
                    F_C[s].setZero();
                    K_CEM[s].setZero();
                    K_CE[s].setZero();
                    K_EMC[s].setZero();
                    for (size_t s2 = 0; s2 < nSpecies; s2++){
                        K_CC[s][s2].setZero();
                    }
                }
                WLumped.setZero();

                for (size_t ip = 0; ip < ipcount; ip++){
                    // no integration here, all done via lumped integration
                    WLumped += w[ip]*Nc[ip];
                }

                //lumped integration
                for (size_t n = 0; n < nNodes_C; n++){  
                    for (size_t s = 0; s < nSpecies; s++){
                        if (ChemPotBased){
                            cLoc[s] = std::exp(C[s](n));
                            dC_dDof[s](n) = std::max(1.0e-4,std::exp(C[s](n)));
                        } else {
                            cLoc[s] = C[s](n);
                            dC_dDof[s](n) = 1.0;
                        }
                        //std::cout << s << " " << cLoc[s] << "   "<< dC_dDof[s](n) << "\n";
                    }


                    for (size_t i = 0; i < nSpecies; i++){
                            cLoc[i] = C[i](n);
                    }  
                    eLoc = E(n);
                    for (size_t r = 0; r < nReactions; r++){
                        for (size_t surf = 0; surf < nReactionTypes; surf++){
                            hasReaction = Reactions[r].GetRates(SurfaceTypes[surf], cLoc, eLoc, EM, ReactionFlux, dReactionFlux_dC, dReactionFlux_dE, iReact, diReact_dE, diReact_dC);
                            if (hasReaction){
                                for (size_t s = 0; s < nSpecies; s++){
                                    F_C[s](n) += -WLumped(n)*ReactionFlux[s];
                                    K_CE[s](n,n) += -WLumped(n)*dReactionFlux_dE[s]*-1.0;
                                    for (size_t s2 = 0; s2 < nSpecies; s2++){
                                        K_CC[s][s2](n,n) += -WLumped(n)*dReactionFlux_dC[s][s2]*dC_dDof[s2](n);
                                    }
                                    if (ChargeConservation){
                                        K_CEM[s](n,0) += -WLumped(n)*dReactionFlux_dE[s];
                                    }
                                }

                                if (ChargeConservation){
                                    F_EM += WLumped(n)*iReact*Curr_Scale;
                                    K_EMEM += WLumped(n)*diReact_dE*Curr_Scale;
                                    K_EME(0,n)  += -WLumped(n)*diReact_dE*Curr_Scale;
                                    for (size_t s = 0; s < nSpecies; s++){
                                        K_EMC[s](0,n) += WLumped(n)*diReact_dC[s]*Curr_Scale*dC_dDof[s](n);
                                    }
                                }
                                ReactionCurrents[r] += WLumped(n)*iReact;
                                if (std::abs(iReact)>=ActiveCurrent_Threshold){
                                    ReactionArea[r] += WLumped(n);
                                }
                            }
                        }
                    }
                    if (ChargeConservation){ //some extra damping
                        K_EMEM += WLumped(n)*Curr_Scale*1e-8;
                    }
                }
                //allocation
                for (size_t i = 0; i < nSpecies; i++){
                    VecAdd(f, dofsC[i], F_C[i]);
                    MatAdd(K, dofsC[i], dofsE, K_CE[i]);
                    for (size_t j = 0; j < nSpecies; j++){
                        MatAdd(K, dofsC[i], dofsC[j], K_CC[i][j]);
                    }
                    if (ChargeConservation){
                        MatAdd(K, dofsC[i], dofEM, K_CEM[i]);
                    }
                }
                if (ChargeConservation){
                    VecAdd(f, dofEM, F_EM);
                    MatAdd(K, dofEM, dofEM, K_EMEM);
                    MatAdd(K, dofEM, dofsE, K_EME);
                    for (size_t i = 0; i < nSpecies; i++){
                        MatAdd(K, dofEM, dofsC[i], K_EMC[i]);
                    }
                }
            }
        }
        if (ChargeConservation){
            PetscMPIInt    rank;
            PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
            std::vector<double> i_reduced(nReactions);
            MPI_Allreduce(ReactionCurrents.data(), i_reduced.data(), nReactions, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
            if (rank==0){
                //double EM_Damping = 0.1/physics->timeScheme->dt*Curr_Scale;
                //MatAdd(K,dofEM, dofEM, EM_Damping); 
                double i_err = 0.0;
                double i_norm = 0.0; 
                for (size_t i = 0; i < i_reduced.size(); i++){
                    i_err += i_reduced[i];
                    i_norm += std::abs(i_reduced[i]);
                    //std::cout << i_reduced[i] << "\n";
                }
                
                std::stringstream ss, ss2; ss << "\tEm=" << EM  << " \n";
                Logs.PrintSingle(ss.str(),3);
                ss2 << "\tCurrent mismatch: " << i_err/i_norm << "\n"; 
                Logs.PrintSingle(ss2.str(),3);
            }
        }
    }
}

size_t ElectroSurface::hasTimeData(std::vector<std::string>& DataNames){
    size_t nData = 2*nReactions+1;
    DataNames.resize(nData);

    for (size_t i = 0; i < nReactions; i++){
        DataNames[2*i] = MyName+"/i_"+Reactions[i].ReactionName;
        DataNames[2*i+1] = MyName+"/A_"+Reactions[i].ReactionName;
    }
            
    DataNames[nData-1] = MyName+"/E_m";

    
    return nData;
}

void ElectroSurface::GetTimeData(std::vector<double>& DataValues){
    size_t nData = nReactions;
    DataValues.resize(0);
    std::vector<double> icurr(nReactions), Acurr(nReactions);

    MPI_Allreduce(ReactionCurrents.data(), icurr.data(), nReactions, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    MPI_Allreduce(ReactionArea.data(), Acurr.data(), nReactions, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

    for (size_t i = 0; i < nReactions; i++){
        DataValues.push_back(icurr[i]);
        DataValues.push_back(Acurr[i]);
    }
    
    if (ChargeConservation){
        double EM;
        PetscInt dofEM;
        dofs->getDofForNodes(EMNode, dofType_EM, dofEM);
        EM = physics->StateVectors[Step_C].GetValue(dofEM);
        DataValues.push_back(EM);
    } else {
        DataValues.push_back(E_m);
    }
}

bool ElectroSurface::SaveData(std::string SaveLoc, std::string DataName, size_t ElemGroup, std::vector<std::vector<double>>& Data){

    //save phasefield energy histroy field
    for (size_t area = 0; area < nAreas; area++){
        if (SaveLoc == "Nodes" && ElemGroup == ElemGroupIndices_C[area] && (DataName == "Overpotential")){
            size_t nNodes = mesh->ElementGroups[ElemGroupIndices_C[area]].NNodes_per_elem;
            PetscInt dofEM;
            double EM;
            if (ChargeConservation){
                dofs->getDofForNodes(EMNode, dofType_EM, dofEM);
                EM = physics->StateVectors[Step_C].GetValue(dofEM);
            } else {
                EM = E_m;
            }



            std::vector<PetscInt> dofsE(nNodes);
            std::vector<size_t> El_E(nNodes);
            Eigen::VectorXd E(nNodes);

            std::vector<Eigen::RowVectorXd> NExport(mesh->ElementGroups[ElemGroupIndices_C[area]].BaseElem->NExport.size()); 
            for (size_t i = 0; i < NExport.size(); i++){
                NExport[i].resize(nNodes);
            }
            for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndices_C[area]].NElems; el++){
                    mesh->GetNodesForElem(El_E, ElemGroupIndices_E[area], el);

                    mesh->getExportShape(ElemGroupIndices_E[area], el, El_E, NExport);

                    dofs->getDofForNodes(El_E, dofType_ePot, dofsE);
                    physics->StateVectors[Step_C].GetValues(dofsE, E);

                    for (size_t n = 0; n < mesh->ElementGroups[ElemGroupIndices_C[area]].BaseElem->NExport.size(); n++){
                        double ePot = NExport[n]*E;

                        Data[el][n] = EM - ePot;
                    }

            }
            return true;
        }
        for (size_t r = 0; r < nReactions; r++){
            if (SaveLoc == "Nodes" && ElemGroup == ElemGroupIndices_C[area] && (DataName == Reactions[r].ReactionName)){
                size_t nNodes = mesh->ElementGroups[ElemGroupIndices_C[area]].NNodes_per_elem;

                PetscInt dofEM;
                double EM;
                if (ChargeConservation){
                    dofs->getDofForNodes(EMNode, dofType_EM, dofEM);
                    EM = physics->StateVectors[Step_C].GetValue(dofEM);
                } else {
                    EM = E_m;
                }

                std::vector<double> ReactionFlux(nSpecies);
                std::vector<std::vector<double>> dReactionFlux_dC(nSpecies); for (size_t i = 0; i < nSpecies; i++) dReactionFlux_dC[i].resize(nSpecies);
                std::vector<double> dReactionFlux_dE(nSpecies);
                double diReact_dE; 
                std::vector<double> diReact_dC(nSpecies);
                std::vector<size_t> RelevantSpecies;

                std::vector<PetscInt> dofsE(nNodes);
                std::vector<std::vector<PetscInt>> dofsC(nSpecies); for (size_t s=0; s<nSpecies; s++) dofsC[s].resize(nNodes);

                std::vector<Eigen::RowVectorXd> NExport(mesh->ElementGroups[ElemGroupIndices_C[area]].BaseElem->NExport.size()); 
                for (size_t i = 0; i < NExport.size(); i++){
                    NExport[i].resize(nNodes);
                }
                     
                std::vector<size_t> El_C(nNodes), El_E(nNodes);
                Eigen::VectorXd E(nNodes);
                std::vector<Eigen::VectorXd> C(nSpecies);
                for (size_t i = 0; i < nSpecies; i++){
                    C[i].resize(nNodes);
                }

                for (size_t el = 0; el < mesh->ElementGroups[ElemGroupIndices_C[area]].NElems; el++){
                    mesh->GetNodesForElem(El_C, ElemGroupIndices_C[area], el);
                    mesh->GetNodesForElem(El_E, ElemGroupIndices_E[area], el);

                    mesh->getExportShape(ElemGroupIndices_C[area], el, El_C, NExport);
                    mesh->getExportShape(ElemGroupIndices_E[area], el, El_E, NExport);

                    dofs->getDofForNodes(El_E, dofType_ePot, dofsE);
                    physics->StateVectors[Step_C].GetValues(dofsE, E);

                    for (size_t s = 0; s < nSpecies; s++){
                        dofs->getDofForNodes(El_C, dofTypes_C[s], dofsC[s]);
                        physics->StateVectors[Step_C].GetValues(dofsC[s], C[s]);    
                    }

                    for (size_t n = 0; n < mesh->ElementGroups[ElemGroupIndices_C[area]].BaseElem->NExport.size(); n++){
                        double ePot = NExport[n]*E;
                        std::vector<double> cVec(nSpecies);
                        for (size_t s = 0; s < nSpecies; s++){
                            if (ChemPotBased){
                                cVec[s] = std::exp(NExport[n]*C[s]);
                            } else {
                                cVec[s] = NExport[n]*C[s];
                            }
                        }

                        double current = 0.0;
                        double AddCurrent;
                        for (size_t i = 0; i < InterfaceReactionTypes[area].size(); i++){
                            Reactions[r].GetRates(InterfaceReactionTypes[area][i], cVec, ePot, EM, ReactionFlux, dReactionFlux_dC, dReactionFlux_dE, AddCurrent, diReact_dE, diReact_dC);
                            current = current+AddCurrent;
                        }
                        Data[el][n] = current;
                    }
                }
                return true;
            }
        }
    }
    return false;
}