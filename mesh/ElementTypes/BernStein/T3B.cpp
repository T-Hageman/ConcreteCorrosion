#include "T3B.h"
#include "../../../utility/utility.h"

T3B::T3B(int ipcount1D){
    Name = "T3B";
    NodeCount = 3;
    dim = 2;
    ipcount = getIPSchemeTriangle(&w, &coordsIP, ipcount1D, dim);
    SetUpBaseShapes();
};

T3B::~T3B(){

};

void T3B::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);;

    for (size_t i = 0; i < ipcount; i++){
        Eigen::Vector2d coords = coordsIP[i];
        double s = coords(0), t = coords(1);

        Nbase[i].resize(3);
        Nbase[i](0) = (1-s-t);
        Nbase[i](1) = s;
        Nbase[i](2) = t;

        Gbase[i].resize(2,3);
        Gbase[i](0,0) = -1.0;    Gbase[i](1,0) = -1.0;
        Gbase[i](0,1) = 1.0;              Gbase[i](1,1) = 0.0;
        Gbase[i](0,2) = 0.0;                Gbase[i](1,2) = 1.0;
    }

    NExport.resize(3);
    for (size_t j = 0; j < 3; j++){
        NExport[j].resize(3);
        NExport[j].setZero();
        NExport[j](j) = 1.0;
    }
    PlottingOrder.resize(3); PlottingOrder[0] = 0; PlottingOrder[1] = 1; PlottingOrder[2] = 2;
}

void T3B::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    getShapeGrads_2D(Nout, Gout, wout, coordsNodes);
}