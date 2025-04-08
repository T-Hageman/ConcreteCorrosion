#include "B2_3.h"
#include "../../../utility/utility.h"

B2_3::B2_3(int ipcount1D){
    Name = "B2_3";
    NodeCount = 16;
    requiresData = true;
    init(ipcount1D, 2);
    SetUpBaseShapes();
};

B2_3::~B2_3(){

};

void B2_3::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);
    G2base.resize(ipcount);

    Eigen::VectorXd Nx(4), Ny(4);
    Eigen::VectorXd Gx(4), Gy(4);
    Eigen::VectorXd G2x(4), G2y(4);
    for (size_t i = 0; i < ipcount; i++){
        Eigen::Vector2d coords = coordsIP[i];

        double x = coords(0), y = coords(1);
        BezierBasis(3,x,Nx,Gx,G2x);
        BezierBasis(3,y,Ny,Gy,G2y);

        Eigen::VectorXd Temp; Temp.resize(Nx.size()*Ny.size());
        Nbase[i].resize(NodeCount);
        KronProd(&Temp, &Nx, &Ny);
        Nbase[i] = Temp.transpose();

        Gbase[i].resize(2,NodeCount);
        KronProd(&Temp, &Gx, &Ny);
        Gbase[i](0,Eigen::indexing::all) = Temp;
        KronProd(&Temp, &Nx, &Gy);
        Gbase[i](1,Eigen::indexing::all) = Temp;

        G2base[i].resize(3,NodeCount);
        KronProd(&Temp, &G2x, &Ny);
        G2base[i](0,Eigen::indexing::all) = Temp;
        KronProd(&Temp, &Nx, &G2y);
        G2base[i](1,Eigen::indexing::all) = Temp;
        KronProd(&Temp, &Gx, &Gy);
        G2base[i](2,Eigen::indexing::all) = Temp;
    }

    NExport.resize(4);
    size_t idx = 0;
    for (size_t j = 0; j < 2; j++){
        for (size_t i = 0; i < 2; i++){
            double x = i;
            double y = j;

            BezierBasis(3,x,Nx,Gx,G2x);
            BezierBasis(3,y,Ny,Gy,G2y);

            Eigen::VectorXd Temp; Temp.resize(Nx.size()*Ny.size());
            KronProd(&Temp, &Nx, &Ny);

            NExport[idx].resize(Nx.size()*Ny.size());
            NExport[idx] = Temp;
            idx += 1;
        }
    }
    PlottingOrder.resize(4); PlottingOrder[0] = 0; PlottingOrder[1] = 1; PlottingOrder[2] = 3; PlottingOrder[3] = 2;
}

void B2_3::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData){
    Eigen::Matrix2d J, Jinv;

    for (size_t i = 0; i < ipcount; i++){
        Nout[i] = Nbase[i]*ElemData.transpose();
        Eigen::MatrixXd GB(2,NodeCount); GB = Gbase[i]*ElemData.transpose();

        J = GB*coordsNodes;
        J.transposeInPlace();
        Jinv = J.inverse();
        wout[i] = w[i]*std::abs(J.determinant());

        for (size_t j = 0; j < NodeCount; j++){
            Gout[i](Eigen::indexing::all,j) = Jinv*GB(Eigen::indexing::all,j);
        }
    }
}

void B2_3::getShapeGrads2(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData){
    Eigen::Matrix2d J, Jinv;
    Eigen::Matrix3d J2, J2inv;
    Eigen::Matrix<double, 3, 2> J2Part2;

    for (size_t i = 0; i < ipcount; i++){
        Nout[i] = Nbase[i]*ElemData.transpose();

        Eigen::MatrixXd GB(2,NodeCount); GB = Gbase[i]*ElemData.transpose();
        Eigen::MatrixXd G2B(3,NodeCount); G2B = G2base[i]*ElemData.transpose();

        J = GB*coordsNodes;
        J2Part2 = G2B*coordsNodes;

        J2(0,0) = J(0,0)*J(0,0);
        J2(0,1) = J(0,1)*J(0,1);
        J2(0,2) = 2*J(0,0)*J(0,1);

        J2(1,0) = J(1,0)*J(1,0);
        J2(1,1) = J(1,1)*J(1,1);
        J2(1,2) = 2*J(1,0)*J(1,1);

        J2(2,0) = J(0,0)*J(1,0);
        J2(2,1) = J(0,1)*J(1,1);
        J2(2,2) = J(1,0)*J(0,1) + J(1,1)*J(0,0);

        Jinv = J.inverse();
        J2inv= J2.inverse();

        wout[i] = w[i]*std::abs(J.determinant());
        for (size_t j = 0; j < NodeCount; j++){
            Gout[i](Eigen::indexing::all,j) = Jinv*GB(Eigen::indexing::all,j);
            G2out[i](Eigen::indexing::all,j)= J2inv*(G2B(Eigen::indexing::all,j) - J2Part2*GB(Eigen::indexing::all,j));
        }        
    }
}

void B2_3::getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData){
    for (size_t i = 0; i < ipcount; i++){
        coordsIP(i,Eigen::indexing::all) = (Nbase[i]*ElemData.transpose())*coordsNodes;
    }
}