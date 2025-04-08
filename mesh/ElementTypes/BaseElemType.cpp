#include "BaseElemType.h"

BaseElemType::BaseElemType(){
    Name = "BaseElem";
    NodeCount = 0;
    requiresData = false;
    requiresNodeData = false;
};

BaseElemType::~BaseElemType(){

};

/// @brief Initializes element, setting up parametric shape function
/// @param ipcount_in Number of integration points per physical dimension
/// @param dim_in Dimension this object exists in
void BaseElemType::init(int ipcount_in, int dim_in){
    dim = dim_in;
    ipcount = getIPScheme(&w, &coordsIP, ipcount_in, dim);

};

/// @brief Evaluate shape functions, gradients and integration weights for an element described by the nodal coordinates coordsNodes 
/// @param Nout output: Element shape functions, N[ip](1:nNodes)
/// @param Gout output: Element gradients, G[ip](x/y, 1:nNodes)
/// @param wout output: Integration weights, w[ip]
/// @param coordsNodes input: Nodal coordinates (1:nNodes, x/y)
void BaseElemType::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    throw std::invalid_argument("Element type does not allow for evaluating shape functions: " + Name);
}

/// @brief Evaluate shape functions, gradients and integration weights for an element described by the nodal coordinates coordsNodes 
/// @param Nout output: Element shape functions, N[ip](1:nNodes)
/// @param Gout output: Element gradients, G[ip](x/y, 1:nNodes)
/// @param wout output: Integration weights, w[ip]
/// @param coordsNodes input: Nodal coordinates (1:nNodes, x/y)
/// @param ElemData input: Element data associated with this element
void BaseElemType::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData){
    throw std::invalid_argument("Element type does not allow for evaluating shape functions: " + Name + "through providing element data\n");
}

/// @brief Evaluate shape functions, gradients and integration weights for an element described by the nodal coordinates coordsNodes 
/// @param Nout output: Element shape functions, N[ip](1:nNodes)
/// @param Gout output: Element gradients, G[ip](x/y, 1:nNodes)
/// @param wout output: Integration weights, w[ip]
/// @param coordsNodes input: Nodal coordinates (1:nNodes, x/y)
/// @param ElemData input: Element data associated with this element
/// @param NodeData input: Node data associated with this element
void BaseElemType::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    throw std::invalid_argument("Element type does not allow for evaluating shape functions: " + Name + "through providing element and node data\n");
}

/// @brief Evaluate shape functions, gradients and integration weights for an element described by the nodal coordinates coordsNodes 
/// @param Nout output: Element shape functions, N[ip](1:nNodes)
/// @param Gout output: Element gradients, G[ip](x/y, 1:nNodes)
/// @param G2out output: Element second gradients, G[ip](xx/yy/xy, 1:nNodes)
/// @param wout output: Integration weights, w[ip]
/// @param coordsNodes input: Nodal coordinates (1:nNodes, x/y)
void BaseElemType::getShapeGrads2(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    throw std::invalid_argument("Element type does not allow for evaluating 2nd order gradient functions: " + Name + "through providing element data\n");    
}

/// @brief Evaluate shape functions, gradients and integration weights for an element described by the nodal coordinates coordsNodes 
/// @param Nout output: Element shape functions, N[ip](1:nNodes)
/// @param Gout output: Element gradients, G[ip](x/y, 1:nNodes)
/// @param G2out output: Element second gradients, G[ip](xx/yy/xy, 1:nNodes)
/// @param wout output: Integration weights, w[ip]
/// @param coordsNodes input: Nodal coordinates (1:nNodes, x/y)
/// @param ElemData input: Element data associated with this element
void BaseElemType::getShapeGrads2(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData){
    throw std::invalid_argument("Element type does not allow for evaluating 2nd order gradient functions: " + Name + "through providing element data\n");    
}

/// @brief Evaluate shape functions, gradients and integration weights for an element described by the nodal coordinates coordsNodes 
/// @param Nout output: Element shape functions, N[ip](1:nNodes)
/// @param Gout output: Element gradients, G[ip](x/y, 1:nNodes)
/// @param G2out output: Element second gradients, G[ip](xx/yy/xy, 1:nNodes)
/// @param wout output: Integration weights, w[ip]
/// @param coordsNodes input: Nodal coordinates (1:nNodes, x/y)
/// @param ElemData input: Element data associated with this element
/// @param NodeData input: Node data associated with this element
void BaseElemType::getShapeGrads2(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    throw std::invalid_argument("Element type does not allow for evaluating 2nd order gradient functions: " + Name + "through providing element and node data\n");    
}


/// @brief Obtain coordinates of integration points
/// @param coordsIP output: Coordinates of integration points 
/// @param coordsNodes inputs: Nodal coordinates
void BaseElemType::getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes){ 
    for (size_t i = 0; i < ipcount; i++){
        coordsIP(i,Eigen::indexing::all) = Nbase[i]*coordsNodes;
    }
}

/// @brief Obtain coordinates of integration points
/// @param coordsIP output: Coordinates of integration points 
/// @param coordsNodes inputs: Nodal coordinates
/// @param ElemData input: Element data associated with this element
void BaseElemType::getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData){
    throw std::invalid_argument("Element type does not allow for evaluating ip coordinates: " + Name + "through providing element data\n");    
}

/// @brief Obtain coordinates of integration points
/// @param coordsIP output: Coordinates of integration points 
/// @param coordsNodes inputs: Nodal coordinates
/// @param ElemData input: Element data associated with this element
/// @param NodeData input: Node data associated with this element
void BaseElemType::getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    throw std::invalid_argument("Element type does not allow for evaluating ip coordinates: " + Name + "through providing element and node data\n");    
}

void BaseElemType::getNormals(std::vector<Eigen::VectorXd> &normals, Eigen::MatrixXd &coordsNodes){
    throw std::invalid_argument("Element type does not allow for evaluating normals: " + Name + "\n");  
}
void BaseElemType::getNormals(std::vector<Eigen::VectorXd> &normals, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData){
    throw std::invalid_argument("Element type does not allow for evaluating normals: " + Name + "through providing element data\n");  
}
void BaseElemType::getNormals(std::vector<Eigen::VectorXd> &normals, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    throw std::invalid_argument("Element type does not allow for evaluating normals: " + Name + "through providing element and node data\n");  
}

/// @brief Obtain values of shape functions at pre-determined points (usually nodes) for exporting results
/// @param Nout Value of shape functions at export points
void BaseElemType::getExportShape(std::vector<Eigen::RowVectorXd> &Nout){
    for (size_t i = 0; i < NExport.size(); i++){
        Nout[i] = NExport[i];
    }
}

/// @brief Obtain values of shape functions at pre-determined points (usually nodes) for exporting results
/// @param Nout Value of shape functions at export points
/// @param ElemData input: Element data associated with this element
void BaseElemType::getExportShape(std::vector<Eigen::RowVectorXd> &Nout, Eigen::MatrixXd &ElemData){
    for (size_t i = 0; i < NExport.size(); i++){
        Nout[i] = NExport[i]*ElemData.transpose();
    }
}

/// @brief Obtain values of shape functions at pre-determined points (usually nodes) for exporting results
/// @param Nout Value of shape functions at export points
/// @param ElemData input: Element data associated with this element
/// @param NodeData input: Node data associated with this element
void BaseElemType::getExportShape(std::vector<Eigen::RowVectorXd> &Nout, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    throw std::invalid_argument("Element type does not allow for evaluating export shape: " + Name + "through providing element and node data\n");   
}

/// @brief Helper function to perform parametric transform for 2D plane shape functions
/// @param Nout 
/// @param Gout 
/// @param wout 
/// @param coordsNodes 
void BaseElemType::getShapeGrads_2D(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    Eigen::Matrix2d J, Jinv;

    for (size_t i = 0; i < ipcount; i++){
        Nout[i] = Nbase[i];
        J = Gbase[i]*coordsNodes;
        J.transposeInPlace();
        Jinv = J.inverse();
        wout[i] = w[i]*std::abs(J.determinant());

        for (size_t j = 0; j < NodeCount; j++){
            Gout[i](Eigen::indexing::all,j) = Jinv*Gbase[i](Eigen::indexing::all,j);
        }
    }
}

/// @brief Helper function to perform parametric transform for 2D line shape functions
/// @param Nout 
/// @param Gout 
/// @param wout 
/// @param coordsNodes 
void BaseElemType::getShapeGrads_2D_Line(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    Eigen::RowVector2d J, tang, normal;
    double JDir;

    for (size_t i = 0; i < ipcount; i++){
        Nout[i] = Nbase[i];
        J = Gbase[i]*coordsNodes;
        tang(0) = J(0); tang(1) = J(1); tang /= std::sqrt(tang.dot(tang));
        normal(0) = tang(1); normal(1) = -tang(0);
        JDir = tang.dot(J);
        wout[i] = w[i]*abs(JDir);
        for (size_t j = 0; j < NodeCount; j++){
            Gout[i](Eigen::indexing::all,j) = (1/JDir)*Gbase[i](Eigen::indexing::all,j);
        }
    }
}

/// @brief Creates Bernstein basic functions
/// @param p order of shape function
/// @param x_ip point to evaluate shape function on
/// @param N output: Shape function values in x_ip
/// @param G output: Shape function gradients in x_ip
/// @param G2 output: Shape function second gradients in x_ip
void BezierBasis(int p, double x_ip, Eigen::VectorXd& N, Eigen::VectorXd& G, Eigen::VectorXd& G2){
    N.resize(p+1);
    G.resize(p+1);
    G2.resize(p+1);

    for (int i = 0; i < p+1; i++){
        double preFactor = fact(p)/fact(p-i)/fact(i);
        N(i) = preFactor*std::pow(1.0-x_ip,std::max(0,p-i))*std::pow(x_ip, std::max(0,i));

        G(i) = preFactor*( (i-p)*std::pow(1.0-x_ip,std::max(0,int(p)-i-1))*std::pow(x_ip, std::max(0,i))
                           +i*std::pow(1.0-x_ip,p-i)*std::pow(x_ip, std::max(0,i-1)));

        G2(i) = preFactor*(2*i      *std::pow(x_ip,std::max(0,i-1))*(i-p)*std::pow(1.0-x_ip,std::max(0,p-i-1))
                           + i      *std::pow(x_ip,std::max(0,i-2))*(i-1)*std::pow(1.0-x_ip,std::max(0,p - i)) 
                           + (i-p+1)*std::pow(x_ip,std::max(0,i))  *(i-p)*std::pow(1.0-x_ip,std::max(0,p-i-2)) );
    }
}

/// @brief Factorial, n!
/// @param n 
/// @return 
int fact(int n){

     return (n==0) || (n==1) ? 1 : n* fact(n-1);
}
