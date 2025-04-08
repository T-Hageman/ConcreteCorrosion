#include "B1_3.h"
#include "../../../utility/utility.h"

B1_3::B1_3(int ipcount1D){
    Name = "B1_3";
    NodeCount = 4;
    requiresData = true;
    init(ipcount1D, 1);
    SetUpBaseShapes();
};

B1_3::~B1_3(){

};

void B1_3::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);
    G2base.resize(ipcount);

    Eigen::VectorXd Nx(4);
    Eigen::VectorXd Gx(4);
    Eigen::VectorXd G2x(4);
    for (size_t i = 0; i < ipcount; i++){
        double x = coordsIP[i](0);

        BezierBasis(3,x,Nx,Gx,G2x);

        Nbase[i].resize(4);
        Nbase[i] = Nx;

        Gbase[i].resize(1,4);
        Gbase[i](0,Eigen::indexing::all) = Gx;

        G2base[i].resize(1,4);
        G2base[i](0,Eigen::indexing::all) = G2x;
    }

    NExport.resize(2);
    size_t idx = 0;
    for (size_t i = 0; i < 2; i++){
        double x = i;
            
        NExport[idx].resize(Nx.size());
        BezierBasis(3,x,Nx,Gx,G2x);
        NExport[idx] = Nx;
        idx += 1;
    }
    PlottingOrder.resize(2); PlottingOrder[0] = 0; PlottingOrder[1] = 1; 
}

void B1_3::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData){
    Eigen::RowVector2d J, tang, normal;
    double JDir;

    for (size_t i = 0; i < ipcount; i++){
        Nout[i] = Nbase[i]*ElemData.transpose();

        Eigen::MatrixXd GB(1,4); GB = Gbase[i]*ElemData.transpose();

        J = GB*coordsNodes;
        tang(0) = J(0); tang(1) = J(1); tang /= std::sqrt(tang.dot(tang));
        normal(0) = tang(1); normal(1) = -tang(0);
        JDir = tang.dot(J);
        wout[i] = w[i]*abs(JDir);
        for (size_t j = 0; j < NodeCount; j++){
            Gout[i](Eigen::indexing::all,j) = (1/JDir)*GB(Eigen::indexing::all,j);
        }
    }
}

void B1_3::getNormals(std::vector<Eigen::VectorXd> &normals, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData){
    Eigen::RowVector2d J, tang, normal;
    Eigen::Vector2d CoordsIp;
    Eigen::RowVectorXd N(NodeCount);
    double JDir;

    for (size_t i = 0; i < ipcount; i++){
        N = Nbase[i]*ElemData.transpose();
        CoordsIp = N*coordsNodes;
        Eigen::MatrixXd GB(1,3); GB = Gbase[i]*ElemData.transpose();

        J = GB*coordsNodes;
        tang(0) = J(0); tang(1) = J(1); tang /= std::sqrt(tang.dot(tang));
        normal(0) = tang(1); normal(1) = -tang(0);
        normals[i] = normal.normalized();

        //ensure outwards Normal
        bool Flip = false;

        double distFromZero = CoordsIp.norm();
        CoordsIp += normals[i]*distFromZero*1e-3;
        double distFromZero2 = CoordsIp.norm();

        if (distFromZero>distFromZero2){
            Flip = true;
        }

        if (Flip){
            normals[i] = -normals[i];
        }
    }
}

void B1_3::getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData){
    for (size_t i = 0; i < ipcount; i++){
        coordsIP(i,Eigen::indexing::all) = Nbase[i]*ElemData.transpose()*coordsNodes;
    }
}