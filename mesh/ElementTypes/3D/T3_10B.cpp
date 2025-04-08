#include "T3_10B.h"
#include "../../../utility/utility.h"

T3_10B::T3_10B(int ipcount1D){
    Name = "T3_10B";
    NodeCount = 10;
    dim = 3;
    ipcount = getIPSchemeTriangle(&w, &coordsIP, ipcount1D, dim);
    SetUpBaseShapes();
};

T3_10B::~T3_10B(){

};

void T3_10B::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);
    G2base.resize(ipcount);

    for (size_t i = 0; i < ipcount; i++){
        Eigen::Vector3d coords = coordsIP[i];

        double x = coords(0), y = coords(1), z = coords(2), v = 1.0-x-y-z;

        Nbase[i].resize(NodeCount);
        Gbase[i].resize(3,NodeCount);
        G2base[i].resize(6,NodeCount);

        Nbase[i](0) = v*v;
        Nbase[i](1) = x*x;
        Nbase[i](2) = y*y;
        Nbase[i](3) = z*z;
        Nbase[i](4) = 2*v*x;
        Nbase[i](5) = 2*x*y;
        Nbase[i](6) = 2*y*v;
        Nbase[i](7) = 2*v*z;
        Nbase[i](8) = 2*x*z;
        Nbase[i](9) = 2*y*z;

        Gbase[i](0,0) = -2*v;           Gbase[i](1,0) = -2*v;           Gbase[i](2,0) = -2*v;
        Gbase[i](0,1) = 2*x;            Gbase[i](1,1) = 0;              Gbase[i](2,1) = 0;
        Gbase[i](0,2) = 0;              Gbase[i](1,2) = 2*y;            Gbase[i](2,2) = 0;
        Gbase[i](0,3) = 0;              Gbase[i](1,3) = 0;              Gbase[i](2,3) = 2*z;
        Gbase[i](0,4) = 2*v-2*x;        Gbase[i](1,4) = -2*x;           Gbase[i](2,4) = -2*x;
        Gbase[i](0,5) = 2*y;            Gbase[i](1,5) = 2*x;            Gbase[i](2,5) = 0;
        Gbase[i](0,6) = -2*y;           Gbase[i](1,6) = 2*v-2*y;        Gbase[i](2,6) = -2*y;
        Gbase[i](0,7) = -2*z;           Gbase[i](1,7) = -2*z;           Gbase[i](2,7) = 2*v-2*z;
        Gbase[i](0,8) = 2*z;            Gbase[i](1,8) = 0;              Gbase[i](2,8) = 2*x;
        Gbase[i](0,9) = 0;              Gbase[i](1,9) = 2*z;            Gbase[i](2,9) = 2*y;

        G2base[i] <<    2, 2, 0, 0, -4, 0, 0, 0, 0, 0, 
                        2, 0, 2, 0, 0,  0, -4, 0, 0, 0,
                        2, 0, 0, 2, 0,  0, 0, -4, 0, 0,
                        2, 0, 0, 0, -2, 2, -2, 0, 0, 0, 
                        2, 0, 0, 0, -2, 0, 0, -2, 4, 0,
                        2, 0, 0, 0, 0,  0, -2,-2, 0, 2; 

    }

    NExport.resize(10);
    for (size_t j = 0; j < 10; j++){
        NExport[j].resize(10);
        NExport[j].setZero();
        NExport[j](j) = 1.0;
    }
    PlottingOrder.resize(4); PlottingOrder[0] = 0; PlottingOrder[1] = 1; PlottingOrder[2] = 2; PlottingOrder[3] = 3;
}

void T3_10B::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    Eigen::Matrix3d J, Jinv;

    for (size_t i = 0; i < ipcount; i++){
        Nout[i] = Nbase[i];

        J = Gbase[i]*coordsNodes;
        J.transposeInPlace();
        Jinv = J.inverse();
        //if (Jinv.hasNaN()){
        //    std::cout << "J= " << J << "\n";
        //}
        wout[i] = w[i]*std::abs(J.determinant());

        for (size_t j = 0; j < NodeCount; j++){
            Gout[i](Eigen::indexing::all,j) = Jinv*Gbase[i](Eigen::indexing::all,j);
        }
    }
}
