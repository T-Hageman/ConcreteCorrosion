#include "L2B.h"
#include "../../../utility/utility.h"

L2B::L2B(int ipcount1D){
    Name = "L2B";
    NodeCount = 2;
    init(ipcount1D, 1);
    SetUpBaseShapes();
};

L2B::~L2B(){

};

void L2B::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);;

    for (size_t i = 0; i < ipcount; i++){
        double coords = coordsIP[i](0);
        Eigen::VectorXd Nx; Nx.resize(2);
        Eigen::VectorXd Gx; Gx.resize(2); 

        double x = coords;
        Nx(0) = (1.0-x);  Gx(0) = -1.0;
        Nx(1) = x;  Gx(1) = 1.0;

        Nbase[i] = Nx.transpose();

        Gbase[i].resize(1,3);
        Gbase[i](0,Eigen::indexing::all) = Gx;
    }

    NExport.resize(2);
    NExport[0].resize(2); NExport[0](0) = 1.0; NExport[0](1) = 0.0;
    NExport[1].resize(2); NExport[1](0) = 0.0; NExport[1](1) = 1.0;
    PlottingOrder.resize(2); PlottingOrder[0] = 0; PlottingOrder[1] = 1; 
}

void L2B::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    getShapeGrads_2D_Line(Nout, Gout, wout, coordsNodes);
}