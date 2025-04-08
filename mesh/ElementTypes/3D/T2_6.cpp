#include "T2_6.h"
#include "../../../utility/utility.h"

T2_6::T2_6(int ipcount1D){
    Name = "T2_6B";
    NodeCount = 6;
    dim = 2;
    ipcount = getIPSchemeTriangle(&w, &coordsIP, ipcount1D, 2);
    SetUpBaseShapes();
};

T2_6::~T2_6(){

};

void T2_6::SetUpBaseShapes(){
    Nbase.resize(ipcount);
    Gbase.resize(ipcount);;

    for (size_t i = 0; i < ipcount; i++){
        Eigen::Vector2d coords = coordsIP[i];
        double s = coords(0), t = coords(1);

        Nbase[i].resize(6);
        Nbase[i](0) = (1-s-t)*(1.0-2*s-2*t);
        Nbase[i](1) = s*(2*s-1.0);
        Nbase[i](2) = t*(2*t-1.0);
        Nbase[i](3) = 4.0*s*(1.0-s-t);
        Nbase[i](4) = 4.0*s*t;
        Nbase[i](5) = 4.0*t*(1.0-s-t);

        Gbase[i].resize(2,6);
        Gbase[i](0,0) = 4.0*s+4.0*t-3.0;    Gbase[i](1,0) = 4.0*s+4.0*t-3.0;
        Gbase[i](0,1) = 4.0*s-1.0;          Gbase[i](1,1) = 0.0;
        Gbase[i](0,2) = 0.0;                Gbase[i](1,2) = 4.0*t-1.0;
        Gbase[i](0,3) = 4.0-4.0*t-8.0*s;    Gbase[i](1,3) = -4.0*s;
        Gbase[i](0,4) = 4.0*t;              Gbase[i](1,4) = 4.0*s;
        Gbase[i](0,5) = -4.0*t;             Gbase[i](1,5) = 4.0-8.0*t-4.0*s;
    }

    NExport.resize(6);
    for (size_t j = 0; j < 6; j++){
        NExport[j].resize(6);
        NExport[j].setZero();
        NExport[j](j) = 1.0;
    }
    PlottingOrder.resize(3); PlottingOrder[0] = 0; PlottingOrder[1] = 1; PlottingOrder[2] = 2;
}

void T2_6::getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes){
    Eigen::RowVector3d tang, tang2;
    Eigen::Matrix<double,2,3> J;
    Eigen::Vector3d JDir;

    for (size_t i = 0; i < ipcount; i++){
        Nout[i] = Nbase[i];

        Eigen::MatrixXd GB(2,6); GB = Gbase[i];

        J = GB*coordsNodes;

        tang(0)  = J(0,0); tang(1)  = J(0,1); tang(2)  = J(0,2); 
        tang2(0) = J(1,0); tang2(1) = J(1,1); tang2(2) = J(1,2); 

        JDir = tang.cross(tang2);

        wout[i] = w[i]*abs(JDir.norm());
        Gout[i].setZero(); //not yet implemented
        //for (size_t j = 0; j < NodeCount; j++){
        //    Gout[i](Eigen::all,j) = (1/JDir)*GB(Eigen::all,j);
        //}
    }
}

