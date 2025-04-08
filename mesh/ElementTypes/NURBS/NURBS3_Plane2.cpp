#include "NURBS3_Plane2.h"
#include "../../../utility/utility.h"

NURBS3_Plane2::NURBS3_Plane2(int ipcount1D){
    Name = "NURBS3_Plane2";
    NodeCount = 9;
    requiresData = true;
    requiresNodeData = true;
    init(ipcount1D, 2);
    SetUpBaseShapes();
};

NURBS3_Plane2::~NURBS3_Plane2(){

};

void NURBS3_Plane2::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);
    G2base.resize(ipcount);

    Eigen::VectorXd Nx(3), Ny(3);
    Eigen::VectorXd Gx(3), Gy(3);
    Eigen::VectorXd G2x(3), G2y(3);
    for (size_t i = 0; i < ipcount; i++){
        Eigen::Vector2d coords = coordsIP[i];

        double x = coords(0), y = coords(1);
        BezierBasis(2,x,Nx,Gx,G2x);
        BezierBasis(2,y,Ny,Gy,G2y);

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

            BezierBasis(2,x,Nx,Gx,G2x);
            BezierBasis(2,y,Ny,Gy,G2y);

            Eigen::VectorXd Temp; Temp.resize(Nx.size()*Ny.size());
            KronProd(&Temp, &Nx, &Ny);

            NExport[idx].resize(Nx.size()*Ny.size());
            NExport[idx] = Temp;
            idx += 1;
        }
    }
    PlottingOrder.resize(4); PlottingOrder[0] = 0; PlottingOrder[1] = 1; PlottingOrder[2] = 3; PlottingOrder[3] = 2;
}

void NURBS3_Plane2::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    Eigen::MatrixXd J(2,3);
    Eigen::RowVector3d tang1, tang2, normal;
    double JDir;

    Eigen::DiagonalMatrix<double, 9> W(NodeData);
    Eigen::VectorXd CB(NodeCount);
    Eigen::MatrixXd GB(NodeCount,2), GXi(NodeCount,2); 
    Eigen::Vector2d dW_dXi; 
    double wFun;

    for (size_t i = 0; i < ipcount; i++){
        CB = ElemData*Nbase[i].transpose();
        wFun = CB.transpose()*NodeData;

        Nout[i] = (W*CB/wFun).transpose();

        GB = ElemData*Gbase[i].transpose();
        dW_dXi = GB.transpose()*NodeData;
        GXi = W*(GB/wFun - CB*dW_dXi.transpose()/wFun/wFun);

        J = GXi.transpose()*coordsNodes;
        tang1 = J(0,Eigen::indexing::all); //tang1.normalize();
        tang2 = J(1,Eigen::indexing::all); //tang2.normalize();
        normal = tang1.cross(tang2);
 
        wout[i] = w[i]*abs(normal.norm());

        // GRADIENT NOT IMPLEMENTED
        //for (size_t j = 0; j < NodeCount; j++){
        //    Gout[i](Eigen::indexing::all,j) = (1/JDir)*GB(Eigen::indexing::all,j);
        //}
    }
}

void NURBS3_Plane2::getNormals(std::vector<Eigen::VectorXd> &normals, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    Eigen::MatrixXd J(2,3);
    Eigen::Vector3d tang1, tang2;
    double JDir;

    Eigen::DiagonalMatrix<double, 9> W(NodeData);
    Eigen::VectorXd CB(NodeCount);
    Eigen::MatrixXd GB(NodeCount,2), GXi(NodeCount,2); 
    Eigen::Vector2d dW_dXi; 
    double wFun;
    Eigen::RowVectorXd N(NodeCount);

    Eigen::Vector3d CoordsIp;

    for (size_t i = 0; i < ipcount; i++){
        CB = ElemData*Nbase[i].transpose();
        wFun = CB.transpose()*NodeData;
        N = (W*CB/wFun).transpose();

        CoordsIp = N*coordsNodes;

        GB = ElemData*Gbase[i].transpose();
        dW_dXi = GB.transpose()*NodeData;
        GXi = W*(GB/wFun - CB*dW_dXi.transpose()/wFun/wFun);

        J = GXi.transpose()*coordsNodes;
        tang1 = J(0,Eigen::indexing::all);
        tang2 = J(1,Eigen::indexing::all);
        normals[i] = tang1.cross(tang2); normals[i].normalize();
        
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

void NURBS3_Plane2::getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    for (size_t i = 0; i < ipcount; i++){
        coordsIP(i,Eigen::indexing::all) = (Nbase[i]*ElemData.transpose())*coordsNodes;
    }
}

void NURBS3_Plane2::getExportShape(std::vector<Eigen::RowVectorXd> &Nout, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData){
    for (size_t i = 0; i < NExport.size(); i++){
        Nout[i] = NExport[i]*ElemData.transpose();
    }
}