#include "T3_10.h"
#include "../../../utility/utility.h"

T3_10::T3_10(int ipcount1D){
    Name = "T3_10";
    NodeCount = 10;
    dim = 3;
    ipcount = getIPSchemeTriangle(&w, &coordsIP, ipcount1D, dim);
    SetUpBaseShapes();
};

T3_10::~T3_10(){

};

void T3_10::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);
    G2base.resize(ipcount);

    for (size_t i = 0; i < ipcount; i++){
        Eigen::Vector3d coords = coordsIP[i];

        double x = coords(0), y = coords(1), z = coords(2);

        Nbase[i].resize(NodeCount);
        Gbase[i].resize(3,NodeCount);
        G2base[i].resize(6,NodeCount);

        Nbase[i](0) = (x+y+z-1)*(2*x+2*y+2*z-1);
        Nbase[i](1) = x*(2*x-1);
        Nbase[i](2) = y*(2*y-1);
        Nbase[i](3) = z*(2*z-1);
        Nbase[i](4) = -x*(4*x+4*y+4*z-4);
        Nbase[i](5) = 4*x*y;
        Nbase[i](6) = -4*y*(x+y+z-1);
        Nbase[i](7) = -z*(4*x+4*y+4*z-4);
        Nbase[i](8) = 4*x*z;
        Nbase[i](9) = 4*y*z;

        Gbase[i](0,0) = 4*x+4*y+4*z-3;  Gbase[i](1,0) = 4*x+4*y+4*z-3;  Gbase[i](2,0) = 4*x+4*y+4*z-3;
        Gbase[i](0,1) = 4*x-1;          Gbase[i](1,1) = 0;              Gbase[i](2,1) = 0;
        Gbase[i](0,2) = 0;              Gbase[i](1,2) = 4*y-1;          Gbase[i](2,2) = 0;
        Gbase[i](0,3) = 0;              Gbase[i](1,3) = 0;              Gbase[i](2,3) = 4*z-1;
        Gbase[i](0,4) = 4-4*y-4*z-8*x;  Gbase[i](1,4) = -4*x;           Gbase[i](2,4) = -4*x;
        Gbase[i](0,5) = 4*y;            Gbase[i](1,5) = 4*x;            Gbase[i](2,5) = 0;
        Gbase[i](0,6) = -4*y;           Gbase[i](1,6) = 4-8*y-4*z-4*x;  Gbase[i](2,6) = -4*y;
        Gbase[i](0,7) = -4*z;           Gbase[i](1,7) = -4*z;           Gbase[i](2,7) = 4-4*y-8*z-4*x;
        Gbase[i](0,8) = 4*z;            Gbase[i](1,8) = 0;              Gbase[i](2,8) = 4*x;
        Gbase[i](0,9) = 0;              Gbase[i](1,9) = 4*z;            Gbase[i](2,9) = 4*y;

        G2base[i] <<    4, 4, 0, 0, -8, 0, 0, 0, 0, 0, 
                        4, 0, 4, 0, 0,  0, -8, 0, 0, 0,
                        4, 0, 0, 4, 0,  0, 0, -8, 0, 0,
                        4, 0, 0, 0, -4, 4, -4, 0, 0, 0, 
                        4, 0, 0, 0, -4, 0, 0, -4, 4, 0,
                        4, 0, 0, 0, 0,  0, -4,-4, 0, 4; 

    }

    NExport.resize(10);
    for (size_t j = 0; j < 10; j++){
        NExport[j].resize(10);
        NExport[j].setZero();
        NExport[j](j) = 1.0;
    }
    PlottingOrder.resize(3); PlottingOrder[0] = 0; PlottingOrder[1] = 1; PlottingOrder[2] = 2; PlottingOrder[3] = 3;
}

void T3_10::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
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
