#include "L3B.h"
#include "../../../utility/utility.h"

L3B::L3B(int ipcount1D){
    Name = "L3B";
    NodeCount = 3;
    init(ipcount1D, 1);
    SetUpBaseShapes();
};

L3B::~L3B(){

};

void L3B::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);;

    for (size_t i = 0; i < ipcount; i++){
        double coords = coordsIP[i](0);
        Eigen::VectorXd Nx; Nx.resize(3);
        Eigen::VectorXd Gx; Gx.resize(3); 

        double x = coords;
        Nx(0) = (1.0-x)*(1.0-x);  Gx(0) = -2.0*(1.0-x);
        Nx(1) = 2*(1.0-x)*x;  Gx(1) = 2*(1.0-x)-2*x;
        Nx(2) = x*x;  Gx(2) = 2*x;

        Nbase[i] = Nx.transpose();

        Gbase[i].resize(1,3);
        Gbase[i](0,Eigen::indexing::all) = Gx;
    }

    NExport.resize(2);
    NExport[0].resize(3); NExport[0](0) = 1.0; NExport[0](1) = 0.0; NExport[0](2) = 0.0;
    NExport[1].resize(3); NExport[1](0) = 0.0; NExport[1](1) = 0.0; NExport[1](2) = 1.0;
    PlottingOrder.resize(2); PlottingOrder[0] = 0; PlottingOrder[1] = 1; 
}

void L3B::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    getShapeGrads_2D_Line(Nout, Gout, wout, coordsNodes);
}