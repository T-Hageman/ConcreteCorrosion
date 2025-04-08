#include "WeakArea.h"
#include "../../Physics/physics.h"

void Register_WeakAreaModel(){
    ModelNames.push_back("WeakArea");
    ModelCreators.push_back(New_WeakAreaModel);
}

BaseModel* New_WeakAreaModel(Physics& My_Physics, std::string MyNameIn){
    return new WeakAreaModel(My_Physics, MyNameIn);
}

WeakAreaModel::WeakAreaModel(Physics& My_Physics, std::string MyName): BaseModel(My_Physics, MyName){
    ModelName = "WeakArea";
};

WeakAreaModel::~WeakAreaModel(){

};

void WeakAreaModel::init(inputData& inputs){
    Setup(inputs);

    std::vector<size_t> UniqueNodes = mesh->ElementGroups[ElementGroupIndex].GetUniqueNodes();
    dofs->AddDofs(UniqueNodes, dofTypes);
}

void WeakAreaModel::load(inputData& inputs, SaveDataFile& data){
    Setup(inputs);
}

void WeakAreaModel::Setup(inputData& inputs){
    std::string EGroupName; inputs.GetRequired(EGroupName, {"Models",MyName,"ElementGroup"});
    ElementGroupIndex = mesh->GetElementGroupIdx(EGroupName);
    
    // Get and Add dofs
    inputs.GetRequired(dofNames, {"Models",MyName,"dofs"});
    dofs->getDofTypesSteps(dofNames, dofTypes, dofSteps);

    //area limits
    inputs.GetRequired(xmin, {"Models",MyName,"xmin"});
    inputs.GetRequired(xmax, {"Models",MyName,"xmax"});

    //other parameters
    inputs.GetRequired(kdummy, {"Models",MyName,"k"});

    inputs.GetRequired(Vals_0,{"Models",MyName,"values"});
    Vals_t.resize(Vals_0.size());
    inputs.GetOptional(Vals_t,{"Models",MyName,"dvalues"});

    LowerTime = -1.0e99;
    UpperTime = 1.0e99;
    inputs.GetOptional(UpperTime, {"Models",MyName,"UpperTime"});

    ExtForce.resize(dofNames.size());

    Lumped = false;
    inputs.GetOptional(Lumped, {"Models",MyName,"Lumped"});
}

void WeakAreaModel::save(SaveDataFile& data){

}

void WeakAreaModel::Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step){
    for (size_t i = 0; i < dofNames.size(); i++){
        if (dofSteps[i] == step){
            ExtForce[i] = 0.0;

            size_t ipcount = mesh->ElementGroups[ElementGroupIndex].BaseElem->ipcount;
            size_t nNodes = mesh->ElementGroups[ElementGroupIndex].NNodes_per_elem;
            Eigen::VectorXd F_el(nNodes);
            std::vector<size_t> Nodes(nNodes);
            Eigen::MatrixXd coordsNodes(nNodes, 2);
            std::vector<Eigen::RowVectorXd> N(ipcount); for (size_t ip=0; ip<ipcount;  ip++) N[ip].resize(nNodes); 
            std::vector<Eigen::MatrixXd> G(ipcount); for (size_t ip=0; ip<ipcount;  ip++) G[ip].resize(1, nNodes); 
            std::vector<double> w(ipcount);
            Eigen::VectorXd X(nNodes);
            Eigen::VectorXd WLumped(nNodes);
            Eigen::MatrixXd K_el(nNodes, nNodes);
            std::vector<PetscInt> dofsX(nNodes);
            Eigen::MatrixXd CoordsIP(ipcount, 2);

            double time = physics->time;
            double dt = physics->timeScheme->dt;
            double Force, targetX;

            if (time>LowerTime && time<UpperTime){
                for (size_t el = 0; el < mesh->ElementGroups[ElementGroupIndex].NElems; el++){
                    mesh->getShapeGrads(ElementGroupIndex, el, Nodes, N, G, w);
                    mesh->GetCoordsForNodes(coordsNodes, Nodes);
                    mesh->getCoordsIP(ElementGroupIndex, el, CoordsIP, coordsNodes);

                    dofs->getDofForNodes(Nodes, dofTypes[i], dofsX);
                    physics->StateVectors[step].GetValues(dofsX, X);

                    F_el.setZero(); 
                    K_el.setZero();
                    WLumped.setZero();
                    for (size_t ip = 0; ip < ipcount; ip++){
                        if (CoordsIP(ip,0)>=xmin && CoordsIP(ip,0)<=xmax){
                            if (Lumped == false){
                                targetX = Vals_0[i]+Vals_t[i]*(time+dt);
                                Force = kdummy * (N[ip]*X-targetX);
                                F_el += -w[ip] * N[ip].transpose()*Force;
                                K_el += -w[ip] * kdummy * N[ip].transpose()*N[ip];

                                ExtForce[i] += w[ip]*Force;
                            } else {
                                WLumped += w[ip] * N[ip];
                            }
                        }
                    }
                    for (size_t n = 0; n < nNodes; n++){
                        targetX = Vals_0[i]+Vals_t[i]*(time+dt);
                        Force = kdummy * (X[n]-targetX);
                        F_el(n)   += -WLumped[n]*Force;
                        K_el(n,n) += -WLumped[n]*kdummy;

                        ExtForce[i] += WLumped[n]*Force;
                    }
                    
                    VecAdd(f, dofsX, F_el);
                    MatAdd(K, dofsX, dofsX, K_el);
                }
            }
        }
    }
}

size_t WeakAreaModel::hasTimeData(std::vector<std::string>& DataNames){
    size_t nData = dofNames.size();
    DataNames.resize(nData);

    for (size_t i = 0; i < nData; i++){
        DataNames[i] = MyName+"/F_"+dofNames[i];
    }
    
    return nData;
}

void WeakAreaModel::GetTimeData(std::vector<double>& DataValues){
    size_t nData = dofNames.size();
    DataValues.resize(nData);
    MPI_Allreduce(ExtForce.data(), DataValues.data(), nData, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
}