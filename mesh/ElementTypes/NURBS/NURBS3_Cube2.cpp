#include "NURBS3_Cube2.h"
#include "../../../utility/utility.h"

NURBS3_Cube2::NURBS3_Cube2(int ipcount1D){
    Name = "NURBS3_Cube2";
    NodeCount = 3*3*3;
    requiresData = true;
    requiresNodeData = true;
    init(ipcount1D, 3);
    SetUpBaseShapes();
};

NURBS3_Cube2::~NURBS3_Cube2(){

};

void NURBS3_Cube2::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);
    G2base.resize(ipcount);

    Eigen::VectorXd Nx(3), Ny(3), Nz(3);
    Eigen::VectorXd Gx(3), Gy(3), Gz(3);
    Eigen::VectorXd G2x(3), G2y(3), G2z(3);
    for (size_t i = 0; i < ipcount; i++){
        Eigen::Vector3d coords = coordsIP[i];

        double x = coords(0), y = coords(1), z = coords(2);
        BezierBasis(2,x,Nx,Gx,G2x);
        BezierBasis(2,y,Ny,Gy,G2y);
        BezierBasis(2,z,Nz,Gz,G2z);

        Eigen::VectorXd Temp; Temp.resize(Nx.size()*Ny.size()*Nz.size());
        Nbase[i].resize(NodeCount);
        KronProd(Temp, Nx, Ny, Nz);
        Nbase[i] = Temp.transpose();

        Gbase[i].resize(3,NodeCount);
        KronProd(Temp, Gx, Ny, Nz);
        Gbase[i](0,Eigen::indexing::all) = Temp;
        KronProd(Temp, Nx, Gy, Nz);
        Gbase[i](1,Eigen::indexing::all) = Temp;
        KronProd(Temp, Nx, Ny, Gz);
        Gbase[i](2,Eigen::indexing::all) = Temp;

        G2base[i].resize(6,NodeCount);
        KronProd(Temp, G2x, Ny, Nz);
        G2base[i](0,Eigen::indexing::all) = Temp;
        KronProd(Temp, Nx, G2y, Nz);
        G2base[i](1,Eigen::indexing::all) = Temp;
        KronProd(Temp, Nx, Ny, G2z);
        G2base[i](2,Eigen::indexing::all) = Temp;
        KronProd(Temp, Gx, Gy, Nz);
        G2base[i](3,Eigen::indexing::all) = Temp;
        KronProd(Temp, Nx, Gy, Gz);
        G2base[i](4,Eigen::indexing::all) = Temp;
        KronProd(Temp, Gx, Ny, Gz);
        G2base[i](5,Eigen::indexing::all) = Temp;

    }

    NExport.resize(8);
    size_t idx = 0;
    for (size_t k = 0; k < 2; k++){
        for (size_t j = 0; j < 2; j++){
            for (size_t i = 0; i < 2; i++){
                double x = i;
                double y = j;
                double z = k;

                BezierBasis(2,x,Nx,Gx,G2x);
                BezierBasis(2,y,Ny,Gy,G2y);
                BezierBasis(2,z,Nz,Gz,G2z);

                Eigen::VectorXd Temp; Temp.resize(Nx.size()*Ny.size()*Nz.size());
                KronProd(Temp, Nx, Ny, Nz);

                NExport[idx].resize(Nx.size()*Ny.size()*Nz.size());
                NExport[idx] = Temp;
                idx += 1;
            }
        }
    }
    PlottingOrder.resize(8); 
    PlottingOrder[0] = 0; PlottingOrder[1] = 1; PlottingOrder[2] = 3; PlottingOrder[3] = 2;
    PlottingOrder[4] = 4; PlottingOrder[5] = 5; PlottingOrder[6] = 7; PlottingOrder[7] = 6;
}

void NURBS3_Cube2::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    Eigen::Matrix3d J, Jinv;
    Eigen::DiagonalMatrix<double, 27> W(NodeData);
    Eigen::VectorXd CB(NodeCount);
    Eigen::MatrixXd GB(NodeCount,3), GXi(NodeCount,3); 
    Eigen::Vector3d dW_dXi; 
    double wFun;

    for (size_t i = 0; i < ipcount; i++){
        CB = ElemData*Nbase[i].transpose();
        wFun = CB.transpose()*NodeData;

        Nout[i] = (W*CB/wFun).transpose();

        GB = ElemData*Gbase[i].transpose();
        dW_dXi = GB.transpose()*NodeData;
        GXi = W*(GB/wFun - CB*dW_dXi.transpose()/wFun/wFun);

        J = GXi.transpose()*coordsNodes;
        Jinv = J.inverse();
        wout[i] = w[i]*std::abs(J.determinant());

        for (size_t j = 0; j < NodeCount; j++){
            Gout[i](Eigen::indexing::all,j) = Jinv*GXi(j,Eigen::indexing::all).transpose();
        }
    }
}

void NURBS3_Cube2::getShapeGrads2(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    Eigen::Matrix3d J, Jinv;
    Eigen::MatrixXd J2(6,6), J2inv(6,6);
    Eigen::Matrix<double, 6, 3> J2Part2;

    Eigen::DiagonalMatrix<double, 27> W(NodeData);
    Eigen::VectorXd CB(NodeCount);
    Eigen::MatrixXd GB(NodeCount,3), GXi(NodeCount,3); 
    Eigen::MatrixXd G2Xi(NodeCount, 6), G2B(NodeCount, 6); 
    Eigen::Vector3d dW_dXi; 
    Eigen::VectorXd dW2_dXi2(6);
    double wFun;

    for (size_t i = 0; i < ipcount; i++){
        CB = ElemData*Nbase[i].transpose();
        wFun = CB.transpose()*NodeData;

        GB = ElemData*Gbase[i].transpose();
        dW_dXi = GB.transpose()*NodeData;
        GXi = W*(GB/wFun - CB*dW_dXi.transpose()/wFun/wFun);

        J = GXi.transpose()*coordsNodes;
        Jinv = J.inverse();

        G2B = ElemData*G2base[i].transpose();
        dW2_dXi2 = G2B.transpose()*NodeData;
        for (size_t ii = 0; ii < 3; ii++){
            for (size_t jj = 0; jj < 3; jj++){
                int ind=-1;
                if (ii==0 && jj==0) ind = 0;
                if (ii==1 && jj==1) ind = 1;
                if (ii==2 && jj==2) ind = 2;

                if (ii==0 && jj==1) ind = 3;
                if (ii==1 && jj==2) ind = 4;
                if (ii==0 && jj==2) ind = 5;

                if (ind>-0.5){
                    G2Xi(Eigen::indexing::all, ind) = W*(
                        G2B(Eigen::indexing::all, ind)/wFun
                        - dW_dXi(jj)*GB(Eigen::indexing::all, ii)/wFun/wFun
                        - dW_dXi(ii)*GB(Eigen::indexing::all, jj)/wFun/wFun
                        - dW2_dXi2(ind)*CB/wFun/wFun
                        + 2*dW_dXi(ii)*dW_dXi(jj)*CB/wFun/wFun/wFun 
                        );
                }
            }
        }
        J2Part2 = G2Xi.transpose()*coordsNodes;

        // https://scicomp.stackexchange.com/questions/25196/implementing-higher-order-derivatives-for-finite-element
        // J(xi/eta/zeta, x/y/z)
        J2(0,0) = J(0,0)*J(0,0);
        J2(0,1) = J(0,1)*J(0,1);
        J2(0,2) = J(0,2)*J(0,2);
        J2(1,0) = J(1,0)*J(1,0);
        J2(1,1) = J(1,1)*J(1,1);
        J2(1,2) = J(1,2)*J(1,2);
        J2(2,0) = J(2,0)*J(2,0);
        J2(2,1) = J(2,1)*J(2,1);
        J2(2,2) = J(2,2)*J(2,2);

        J2(0,3) = 2*J(0,0)*J(0,1);
        J2(0,4) = 2*J(0,1)*J(0,2);
        J2(0,5) = 2*J(0,0)*J(0,2);
        J2(1,3) = 2*J(1,0)*J(1,1);
        J2(1,4) = 2*J(1,1)*J(1,2);
        J2(1,5) = 2*J(1,0)*J(1,2);
        J2(2,3) = 2*J(2,0)*J(2,1);
        J2(2,4) = 2*J(2,1)*J(2,2);
        J2(2,5) = 2*J(2,0)*J(2,2);

        J2(3,0) = J(0,0)*J(1,0);
        J2(3,1) = J(0,1)*J(1,1);
        J2(3,2) = J(0,2)*J(1,2);
        J2(4,0) = J(1,0)*J(2,0);
        J2(4,1) = J(1,1)*J(2,1);
        J2(4,2) = J(1,2)*J(2,2);
        J2(5,0) = J(0,0)*J(2,0);
        J2(5,1) = J(0,1)*J(2,1);
        J2(5,2) = J(0,2)*J(2,2);

        J2(3,3) = J(0,0)*J(1,1) + J(1,0)*J(0,1);
        J2(3,4) = J(0,1)*J(1,2) + J(1,1)*J(0,2);
        J2(3,5) = J(0,0)*J(1,2) + J(1,0)*J(0,2);

        J2(4,3) = J(1,0)*J(2,1) + J(2,0)*J(1,1);
        J2(4,4) = J(1,1)*J(2,2) + J(2,1)*J(1,2);
        J2(4,5) = J(1,0)*J(2,2) + J(2,0)*J(1,2);

        J2(5,3) = J(0,0)*J(2,1) + J(2,0)*J(0,1);
        J2(5,4) = J(0,1)*J(2,2) + J(2,1)*J(0,2);
        J2(5,5) = J(0,0)*J(2,2) + J(2,0)*J(0,2);

        J2inv= J2.inverse();

        wout[i] = w[i]*std::abs(J.determinant());
        Nout[i] = (W*CB/wFun).transpose();
        for (size_t j = 0; j < NodeCount; j++){
            Gout[i](Eigen::indexing::all,j) = Jinv*GXi(j,Eigen::indexing::all).transpose();
            G2out[i](Eigen::indexing::all,j)= J2inv*(G2Xi(j,Eigen::indexing::all).transpose() - J2Part2*GXi(j,Eigen::indexing::all).transpose());
        }        
    }
}

void NURBS3_Cube2::getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    Eigen::VectorXd CB(NodeCount);
    double wFun;
    Eigen::DiagonalMatrix<double, 27> W(NodeData);

    for (size_t i = 0; i < ipcount; i++){
        CB = ElemData*Nbase[i].transpose();
        wFun = CB.transpose()*NodeData;

        coordsIP(i,Eigen::indexing::all) = (W*CB/wFun).transpose()*coordsNodes;
    }
}

void NURBS3_Cube2::getExportShape(std::vector<Eigen::RowVectorXd> &Nout, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    Eigen::VectorXd CB(NodeCount);
    double wFun;
    Eigen::DiagonalMatrix<double, 27> W(NodeData);

    for (size_t i = 0; i < NExport.size(); i++){
        CB = ElemData*NExport[i].transpose();
        wFun = CB.transpose()*NodeData;

        Nout[i] = (W*CB/wFun).transpose();
    }
}