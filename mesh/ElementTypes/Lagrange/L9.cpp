#include "L9.h"
#include "../../../utility/utility.h"

L9::L9(int ipcount1D){
    Name = "L9";
    NodeCount = 9;
    init(ipcount1D, 2);
    SetUpBaseShapes();
};

L9::~L9(){

};

void L9::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);;

    Eigen::VectorXd Nx, Ny; Nx.resize(3); Ny.resize(3);
    Eigen::VectorXd Gx, Gy; Gx.resize(3); Gy.resize(3);
    for (size_t i = 0; i < ipcount; i++){
        Eigen::Vector2d coords = coordsIP[i];

        double x = coords(0), y = coords(1);
        Nx(0) = (x-0.5)*(x-1)/((0-0.5)*(0-1));  Gx(0) = (2*x-1.5)/((0-0.5)*(0-1));
        Nx(1) = (x-0)*(x-1)/((0.5-0)*(0.5-1));  Gx(1) = (2*x-1.0)/((0.5-0)*(0.5-1));
        Nx(2) = (x-0)*(x-0.5)/((1-0)*(1-0.5));  Gx(2) = (2*x-0.5)/((1-0)*(1-0.5));
        Ny(0) = (y-0.5)*(y-1)/((0-0.5)*(0-1));  Gy(0) = (2*y-1.5)/((0-0.5)*(0-1));
        Ny(1) = (y-0)*(y-1)/((0.5-0)*(0.5-1));  Gy(1) = (2*y-1.0)/((0.5-0)*(0.5-1));
        Ny(2) = (y-0)*(y-0.5)/((1-0)*(1-0.5));  Gy(2) = (2*y-0.5)/((1-0)*(1-0.5));

        Eigen::VectorXd Temp; Temp.resize(Nx.size()*Ny.size());

        Nbase[i].resize(9);
        KronProd(&Temp, &Nx, &Ny);
        Nbase[i] = Temp.transpose();

        Gbase[i].resize(2,9);
        KronProd(&Temp, &Gx, &Ny);
        Gbase[i](0,Eigen::indexing::all) = Temp;
        KronProd(&Temp, &Nx, &Gy);
        Gbase[i](1,Eigen::indexing::all) = Temp;
    }

    NExport.resize(4);
    size_t idx = 0;
    for (size_t j = 0; j < 2; j++){
        for (size_t i = 0; i < 2; i++){
            Nx(0) = (1.0-i); Nx(1) = 0.0; Nx(2) = i;
            Ny(0) = (1.0-j); Ny(1) = 0.0; Ny(2) = j;

            Eigen::VectorXd Temp; Temp.resize(Nx.size()*Ny.size());
            KronProd(&Temp, &Nx, &Ny);

            NExport[idx].resize(Nx.size()*Ny.size());
            NExport[idx] = Temp;
            idx += 1;
        }
    }
    PlottingOrder.resize(4); PlottingOrder[0] = 0; PlottingOrder[1] = 1; PlottingOrder[2] = 3; PlottingOrder[3] = 2;
}

void L9::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    getShapeGrads_2D(Nout, Gout, wout, coordsNodes);
}