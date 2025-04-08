#include "L3.h"
#include "../../../utility/utility.h"

L3::L3(int ipcount1D){
    Name = "L3";
    NodeCount = 3;
    init(ipcount1D, 1);
    SetUpBaseShapes();
};

L3::~L3(){

};

void L3::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);;

    for (size_t i = 0; i < ipcount; i++){
        double coords = coordsIP[i](0);
        Eigen::VectorXd Nx; Nx.resize(3);
        Eigen::VectorXd Gx; Gx.resize(3); 

        double x = coords;
        Nx(0) = (x-0.5)*(x-1)/((0-0.5)*(0-1));  Gx(0) = (2*x-1.5)/((0-0.5)*(0-1));
        Nx(1) = (x-0)*(x-1)/((0.5-0)*(0.5-1));  Gx(1) = (2*x-1.0)/((0.5-0)*(0.5-1));
        Nx(2) = (x-0)*(x-0.5)/((1-0)*(1-0.5));  Gx(2) = (2*x-0.5)/((1-0)*(1-0.5));

        Nbase[i] = Nx.transpose();

        Gbase[i].resize(1,3);
        Gbase[i](0,Eigen::indexing::all) = Gx;
    }

    NExport.resize(2);
    NExport[0].resize(3); NExport[0](0) = 1.0; NExport[0](1) = 0.0; NExport[0](2) = 0.0;
    NExport[1].resize(3); NExport[1](0) = 0.0; NExport[1](1) = 0.0; NExport[1](2) = 1.0;
    PlottingOrder.resize(2); PlottingOrder[0] = 0; PlottingOrder[1] = 1;
}

void L3::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    getShapeGrads_2D_Line(Nout, Gout, wout, coordsNodes);
}