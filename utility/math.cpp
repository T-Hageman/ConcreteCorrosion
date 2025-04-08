#include "utility.h"

void ConstructBMat(Eigen::MatrixXd& B, Eigen::MatrixXd& G, int dim){
    size_t nNodes = G.cols();

    B.setZero();
    if (dim==2){ //xx/yy/zz/xy
        B(0, Eigen::seq(0, nNodes - 1)) = G(0, Eigen::indexing::all);
        B(3, Eigen::seq(0, nNodes - 1)) = G(1, Eigen::indexing::all);
        B(1, Eigen::seq(nNodes, 2 * nNodes - 1)) = G(1, Eigen::indexing::all);
        B(3, Eigen::seq(nNodes, 2 * nNodes - 1)) = G(0, Eigen::indexing::all);
    }
    if (dim==3){ //xx/yy/zz/xy/yz/xz
        B(0, Eigen::seq(0, nNodes - 1)) = G(0, Eigen::indexing::all);
        B(3, Eigen::seq(0, nNodes - 1)) = G(1, Eigen::indexing::all);
        B(5, Eigen::seq(0, nNodes - 1)) = G(2, Eigen::indexing::all);

        B(1, Eigen::seq(nNodes, 2 * nNodes - 1)) = G(1, Eigen::indexing::all);
        B(3, Eigen::seq(nNodes, 2 * nNodes - 1)) = G(0, Eigen::indexing::all);
        B(4, Eigen::seq(nNodes, 2 * nNodes - 1)) = G(2, Eigen::indexing::all);

        B(2, Eigen::seq(2*nNodes, 3 * nNodes - 1)) = G(2, Eigen::indexing::all);
        B(4, Eigen::seq(2*nNodes, 3 * nNodes - 1)) = G(1, Eigen::indexing::all);
        B(5, Eigen::seq(2*nNodes, 3 * nNodes - 1)) = G(0, Eigen::indexing::all);
    }
}


/// @brief Kronecker product of two vectors
/// @param result 
/// @param in1 
/// @param in2 
void KronProd(Eigen::VectorXd* result, Eigen::VectorXd* in1, Eigen::VectorXd* in2){
    int sizeInner = in1->size();
    for (long int j = 0; j < in2->size(); j++){
        for (long int i = 0; i < sizeInner; i++){
            result->operator()(i+sizeInner*j) = in1->operator()(i)*in2->operator()(j);
        }   
    }
};

void KronProd(Eigen::VectorXd& result, Eigen::VectorXd& in1, Eigen::VectorXd& in2){
    int sizeInner = in1.size();
    for (long int j = 0; j < in2.size(); j++){
        for (long int i = 0; i < sizeInner; i++){
            result(i+sizeInner*j) = in1(i)*in2(j);
        }   
    }
};

void KronProd(Eigen::VectorXd& result, Eigen::VectorXd& in1, Eigen::VectorXd& in2, Eigen::VectorXd& in3){
    int cnt = 0;
    for (size_t k = 0; k < in3.size(); k++){
        for (size_t j = 0; j < in2.size(); j++){
            for (size_t i = 0; i < in1.size(); i++){
                result(cnt) = in1(i)*in2(j)*in3(k);
                cnt += 1;
            }   
        }
    }
};

double min(std::vector<std::vector<double>> &Matrix){
    if (Matrix.size()==0){
        return 0;
    }
    if (Matrix[0].size()==0){
        return 0;
    }
    
    double minval = Matrix[0][0];
    for (size_t i = 0; i < Matrix.size(); i++){
        for (size_t j = 0; j < Matrix[0].size(); j++){
            if (Matrix[i][j]<minval){
                minval = Matrix[i][j];
            }
        }   
    }
    return minval;
}

double max(std::vector<std::vector<double>> &Matrix){
    if (Matrix.size()==0){
        return 0;
    }
    if (Matrix[0].size()==0){
        return 0;
    }
    
    double maxval = Matrix[0][0];
    for (size_t i = 0; i < Matrix.size(); i++){
        for (size_t j = 0; j < Matrix[0].size(); j++){
            if (Matrix[i][j]>maxval){
                maxval = Matrix[i][j];
            }
        }   
    }
    return maxval;
}

double sgn(double x){
    if (x > 0) return 1;
    if (x < 0) return -1;
    return 0;
}

double sum(std::vector<double>& ToSum){
    double res = std::reduce(ToSum.begin(), ToSum.end());
    return res;
}

size_t sum(std::vector<size_t>& ToSum){
    double res = std::reduce(ToSum.begin(), ToSum.end());
    return res;
}

GhostVector::GhostVector(){

};

GhostVector::~GhostVector(){

}

/// @brief For single core, sets the size of this vector and initializes the mapping
/// @param nGlobal Global size of this vector
void GhostVector::SetSizeSequential(size_t nGlobal){
    PetscCallThrow(VecCreateSeq(PETSC_COMM_WORLD, nGlobal, &DataVector));

    std::vector<size_t> localNums(nGlobal); std::iota(localNums.begin(), localNums.end(),0);

    for (size_t i = 0; i < localNums.size(); i++){
        GlobalToLocal.insert({localNums[i], localNums[i]});
        LocalToGlobal.insert({localNums[i], localNums[i]}); 
    }
}

/// @brief Sets size of this vector, initializes mapping operators, and sets up ghost Nodes
/// @param nLocal Local size of this vector
/// @param nGlobal global size of this vector
/// @param GhostLocs Location of ghost nodes
void GhostVector::SetSize(size_t nLocal, size_t nGlobal, std::vector<PetscInt> &GhostLocs){
    size_t nGhosts = GhostLocs.size();
    std::vector<PetscInt> Ownerrange(2);

    PetscCallThrow(VecCreateGhost(PETSC_COMM_WORLD, nLocal, nGlobal, nGhosts, GhostLocs.data(), &DataVector));
    PetscCallThrow(VecGetOwnershipRange(DataVector, &Ownerrange[0], &Ownerrange[1]));

    std::vector<size_t> globalNums(nLocal+nGhosts); 
    for (size_t i = 0; i < nLocal; i++){
        globalNums[i] = Ownerrange[0]+i;
    }
    Ghosts.resize(nGhosts);
    for (size_t i = 0; i < nGhosts; i++){
        globalNums[nLocal+i] = GhostLocs[i];
        Ghosts[i] = GhostLocs[i];
    }
    std::vector<size_t> localNums(nLocal+nGhosts); std::iota(localNums.begin(), localNums.end(),0);

    for (size_t i = 0; i < localNums.size(); i++){
        GlobalToLocal.insert({globalNums[i], localNums[i]});
        LocalToGlobal.insert({localNums[i], globalNums[i]}); 
    }
}

/// @brief Sets or adds values to ghost vector
/// @param locs global index locations to set
/// @param vals values to set vector to
/// @param operation INSERT_VALUES / ADD_VALUES
void GhostVector::Set(std::vector<PetscInt> &locs, std::vector<double> &vals, InsertMode operation){
    PetscCallThrow(VecSetValues(DataVector, locs.size(), locs.data(), vals.data(), operation));
}

/// @brief Sets or adds a constant value to ghost vector
/// @param locs global index locations to set
/// @param vals value to set all locs of the vector to
/// @param operation INSERT_VALUES / ADD_VALUES
void GhostVector::Set(PetscInt locs, double vals, InsertMode operation){
    PetscCallThrow(VecSetValue(DataVector, locs, vals, operation));
}

/// @brief Zeros ghost vector
void GhostVector::Zero(){
    PetscCallThrow(VecZeroEntries(DataVector));
}

/// @brief start flush of all LOCAL changes to vector
void GhostVector::AssemblyStart(){
    PetscCallThrow(VecAssemblyBegin(DataVector));
}

/// @brief ends flush of all LOCAL changes to vector
void GhostVector::AssemblyEnd(){
    PetscCallThrow(VecAssemblyEnd(DataVector));
}

/// @brief Start sync changes across ghost nodes
/// @param op INSERT_VALUES / ADD_VALUES
void GhostVector::SyncStart(InsertMode op){
    if (op == INSERT_VALUES){
        PetscCallThrow(VecGhostUpdateBegin(DataVector, op, SCATTER_FORWARD));
    } else {
        PetscCallThrow(VecGhostUpdateBegin(DataVector, op, SCATTER_REVERSE));
    }
}

/// @brief End sync changes across ghost nodes
/// @param op INSERT_VALUES / ADD_VALUES
void GhostVector::SyncEnd(InsertMode op){
    if (op == INSERT_VALUES){
        PetscCallThrow(VecGhostUpdateEnd(DataVector, op, SCATTER_FORWARD));
    } else {
        PetscCallThrow(VecGhostUpdateEnd(DataVector, op, SCATTER_REVERSE));
    }
}

/// @brief Obtain values from the ghost vector using global indices
/// @param locs input: global indices to obtain values from
/// @param vals_out output: Values obtained
void GhostVector::GetValues(std::vector<size_t> &locs, Eigen::VectorXd &vals_out){
    Vec LocalVec;
    PetscCallThrow(VecGhostGetLocalForm(DataVector, &LocalVec));

    PetscScalar* V;
    
    PetscCallThrow(VecGetArray(LocalVec, &V));
    for (size_t i = 0; i < locs.size(); i++){
        size_t localNum = GlobalToLocal[locs[i]];
        vals_out[i] = V[localNum];
    }
    PetscCallThrow(VecRestoreArray(LocalVec, &V));
    PetscCallThrow(VecGhostRestoreLocalForm(DataVector,&LocalVec)); 
}

/// @brief Obtain values from the ghost vector using global indices
/// @param locs input: global indices to obtain values from
/// @param vals_out output: Values obtained
void GhostVector::GetValues(std::vector<PetscInt> &locs, Eigen::VectorXd &vals_out){
    Vec LocalVec;
    PetscCallThrow(VecGhostGetLocalForm(DataVector, &LocalVec));

    PetscScalar* V;
    
    PetscCallThrow(VecGetArray(LocalVec, &V));
    for (size_t i = 0; i < locs.size(); i++){
        size_t localNum = GlobalToLocal[locs[i]];
        vals_out[i] = V[localNum];
    }
    PetscCallThrow(VecRestoreArray(LocalVec, &V));
    PetscCallThrow(VecGhostRestoreLocalForm(DataVector,&LocalVec));     
}

/// @brief Obtain values from the ghost vector using global indices
/// @param locs input: global indices to obtain values from
/// @param vals_out output: Values obtained
void GhostVector::GetValues(std::vector<size_t> &locs, std::vector<double> &vals_out){
    Vec LocalVec;
    PetscCallThrow(VecGhostGetLocalForm(DataVector, &LocalVec));

    PetscScalar* V;
    
    PetscCallThrow(VecGetArray(LocalVec, &V));
    for (size_t i = 0; i < locs.size(); i++){
        size_t localNum = GlobalToLocal[locs[i]];
        vals_out[i] = V[localNum];
    }
    PetscCallThrow(VecRestoreArray(LocalVec, &V));
    PetscCallThrow(VecGhostRestoreLocalForm(DataVector,&LocalVec));        
}

/// @brief Gets a single value from the vector
/// @param locs Location in global index
/// @return Value obtained from vector
double GhostVector::GetValue(PetscInt locs){
    Vec LocalVec;
    double val_out;
    PetscCallThrow(VecGhostGetLocalForm(DataVector, &LocalVec));

    PetscScalar* V;
    
    PetscCallThrow(VecGetArray(LocalVec, &V));
    size_t localNum = GlobalToLocal[locs];
    val_out = V[localNum];
    PetscCallThrow(VecRestoreArray(LocalVec, &V));
    PetscCallThrow(VecGhostRestoreLocalForm(DataVector,&LocalVec));       
    return val_out;
}

/// @brief Sets values in vector based on other ghost vector
/// @param ToCopy Ghost vector to copy values from
void GhostVector::Set(GhostVector& ToCopy){
    //PetscCallThrow(VecAYPX(DataVector, 0.0, ToCopy.DataVector));
    PetscCallThrow(VecCopy(ToCopy.DataVector, DataVector));
}

/// @brief Adds contribution to a PETSC matrix object
/// @param K Matrix to add into
/// @param rows Global index rows to add into
/// @param cols Global index columns to add into
/// @param K_el Smaller matrix to add to large matrix
void MatAdd(Mat& K, std::vector<PetscInt> &rows, std::vector<PetscInt> &cols, Eigen::MatrixXd &K_el){
    PetscCallThrow(MatSetValues(K, rows.size(), rows.data(), cols.size(), cols.data(), K_el.data(), ADD_VALUES));
}

void MatAdd(Mat& K, PetscInt &rows, std::vector<PetscInt> &cols, Eigen::MatrixXd &K_el){
    PetscCallThrow(MatSetValues(K, 1, &rows, cols.size(), cols.data(), K_el.data(), ADD_VALUES));
}

void MatAdd(Mat& K, std::vector<PetscInt> &rows, PetscInt &cols, Eigen::MatrixXd &K_el){
    PetscCallThrow(MatSetValues(K, rows.size(), rows.data(), 1, &cols, K_el.data(), ADD_VALUES));
}

void MatAdd(Mat& K, PetscInt &rows, PetscInt &cols, double &K_el){
    PetscCallThrow(MatSetValue(K, rows, cols, K_el, ADD_VALUES));
}

/// @brief Add values to a PETSC vector
/// @param V Vector to add to
/// @param rows Global indices to add to
/// @param F_el Values to add
void VecAdd(Vec& V, std::vector<PetscInt> &rows, Eigen::VectorXd &F_el){
    PetscCallThrow(VecSetValues(V, rows.size(), rows.data(), F_el.data(), ADD_VALUES));
}

/// @brief Add values to a PETSC vector
/// @param V Vector to add to
/// @param rows Global indices to add to
/// @param F_el Values to add
void VecAdd(Vec& V, PetscInt &rows, double &F_el){
    PetscCallThrow(VecSetValue(V, rows, F_el, ADD_VALUES));
}

void GetLog10(std::vector<std::vector<double>> &NodeData){
    for (size_t i = 0; i < NodeData.size(); i++){
        for (size_t j = 0; j < NodeData[0].size(); j++){
            NodeData[i][j] = std::log10(std::max(std::abs(NodeData[i][j]),1.0e-16));
        }
    }
}

size_t Find(std::vector<std::string> array, std::string ToFind, bool& found){
    size_t index;
    auto it = std::find(array.begin(), array.end(), ToFind);
    if (it == array.end()){
        index = 0;
        found = false;
    } else{
        index = std::distance(array.begin(), it);
        found = true;
    }
    return index;
}