#include "L2.h"
#include "../../../utility/utility.h"

L2::L2(int ipcount1D){
    Name = "L2";
    NodeCount = 2;
    init(ipcount1D, 1);
    SetUpBaseShapes();
};

L2::~L2(){

};

void L2::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);;

    for (size_t i = 0; i < ipcount; i++){
        double coords = coordsIP[i](0);
        Eigen::VectorXd Nx; Nx.resize(2);
        Eigen::VectorXd Gx; Gx.resize(2); 

        double x = coords;
        Nx(0) = (x-1.0)/(0.0-1.0);  Gx(0) = (1.0)/(0.0-1.0);
        Nx(1) = (x-0.0)/(1.0-0.0);  Gx(1) = (1.0)/(1.0-0.0);


        Nbase[i] = Nx.transpose();

        Gbase[i].resize(1,2);
        Gbase[i](0,Eigen::indexing::all) = Gx;
    }

    NExport.resize(2);
    NExport[0].resize(2); NExport[0](0) = 1.0; NExport[0](1) = 0.0;
    NExport[1].resize(2); NExport[1](0) = 0.0; NExport[1](1) = 1.0;
    PlottingOrder.resize(2); PlottingOrder[0] = 0; PlottingOrder[1] = 1;
}

void L2::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    getShapeGrads_2D_Line(Nout, Gout, wout, coordsNodes);
}