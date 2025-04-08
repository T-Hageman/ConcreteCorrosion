#include "L4.h"
#include "../../../utility/utility.h"

L4::L4(int ipcount1D){
    Name = "L4";
    NodeCount = 4;
    init(ipcount1D, 2);
    SetUpBaseShapes();
};

L4::~L4(){

};

void L4::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);

    Eigen::VectorXd Nx, Ny; Nx.resize(2); Ny.resize(2);
    Eigen::VectorXd Gx, Gy; Gx.resize(2); Gy.resize(2);
    for (size_t i = 0; i < ipcount; i++){
        Eigen::Vector2d coords = coordsIP[i];

        double x = coords(0), y = coords(1);
        Nx(0) = (x-1.0)/(0.0-1.0);  Gx(0) = (1.0)/(0.0-1.0);
        Nx(1) = (x-0.0)/(1.0-0.0);  Gx(1) = (1.0)/(1.0-0.0);

        Ny(0) = (y-1.0)/(0.0-1.0);  Gy(0) = (1.0)/(0.0-1.0);
        Ny(1) = (y-0.0)/(1.0-0.0);  Gy(1) = (1.0)/(1.0-0.0);

        Eigen::VectorXd Temp; Temp.resize(Nx.size()*Ny.size());

        Nbase[i].resize(4);
        KronProd(&Temp, &Nx, &Ny);
        Nbase[i] = Temp.transpose();

        Gbase[i].resize(2,4);
        KronProd(&Temp, &Gx, &Ny);
        Gbase[i](0,Eigen::indexing::all) = Temp;
        KronProd(&Temp, &Nx, &Gy);
        Gbase[i](1,Eigen::indexing::all) = Temp;
    }

    NExport.resize(4);
    size_t idx = 0;
    for (size_t j = 0; j < 2; j++){
        for (size_t i = 0; i < 2; i++){
            Nx(0) = (1.0-i); Nx(1) = i;
            Ny(0) = (1.0-j); Ny(1) = j;

            Eigen::VectorXd Temp; Temp.resize(Nx.size()*Ny.size());
            KronProd(&Temp, &Nx, &Ny);

            NExport[idx].resize(Nx.size()*Ny.size());
            NExport[idx] = Temp;
            idx += 1;
        }
    }
    PlottingOrder.resize(4); PlottingOrder[0] = 0; PlottingOrder[1] = 1; PlottingOrder[2] = 3; PlottingOrder[3] = 2;
}

void L4::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    getShapeGrads_2D(Nout, Gout, wout, coordsNodes);
}