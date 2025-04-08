#include "isotropicDiff.h"

double GetIsotropicDiffusivity(Eigen::RowVectorXd &Nu,Eigen::VectorXd& uNodes, std::vector<double>& w, double diff){
    double scale;
    double wSum = std::accumulate(w.begin(), w.end(), 0.0);
    scale = std::sqrt(wSum);

    size_t nNodes_Vel = uNodes.size()/2;
    Eigen::MatrixXd U2(nNodes_Vel,2);
    U2(Eigen::indexing::all,0) = uNodes(Eigen::seq(0,nNodes_Vel-1));
    U2(Eigen::indexing::all,1) = uNodes(Eigen::seq(nNodes_Vel,2*nNodes_Vel-1));

    Eigen::Vector2d u;
    u = Nu*U2;

    double uNorm = std::sqrt(u.transpose()*u);

    double Pe_Target = 1.0;

    double diff_art = uNorm*scale/2.0/Pe_Target;

    if (diff_art<=diff){
        return diff;
    } else {
        return diff_art;
    }
}

Eigen::Matrix4d GetIsotropicDiffusivity(Eigen::RowVectorXd &Nu,Eigen::VectorXd& uNodes, std::vector<double>& w, Eigen::Matrix4d& diff){
    Eigen::Matrix4d Diff_Stab; Diff_Stab=diff;

    double scale;
    double wSum = std::accumulate(w.begin(), w.end(), 0.0);
    scale = std::sqrt(wSum);

    size_t nNodes_Vel = uNodes.size()/2;
    Eigen::MatrixXd U2(nNodes_Vel,2);
    U2(Eigen::indexing::all,0) = uNodes(Eigen::seq(0,nNodes_Vel-1));
    U2(Eigen::indexing::all,1) = uNodes(Eigen::seq(nNodes_Vel,2*nNodes_Vel-1));

    Eigen::Vector2d u;
    u = Nu*U2;
    u(0) = abs(u(0));
    u(1) = abs(u(1));

    double Pe_Target = 1.0;

    Eigen::Vector2d diff_art = u*scale/2.0/Pe_Target;

    if (Diff_Stab(0,0)<diff_art(0)){
        Diff_Stab(0,0)=diff_art(0);
    }
    if (Diff_Stab(1,1)<diff_art(1)){
        Diff_Stab(1,1)=diff_art(1);
    }

    return Diff_Stab;
}

Eigen::Matrix4d GetIsotropicDiffusivity(Eigen::RowVectorXd &Nu,Eigen::VectorXd& uNodes, std::vector<double>& w, Eigen::Matrix4d& diff, double rho){
    Eigen::Matrix4d Diff_Stab; Diff_Stab=diff;

    double scale;
    double wSum = std::accumulate(w.begin(), w.end(), 0.0);
    scale = std::sqrt(wSum);

    size_t nNodes_Vel = uNodes.size()/2;
    Eigen::MatrixXd U2(nNodes_Vel,2);
    U2(Eigen::indexing::all,0) = uNodes(Eigen::seq(0,nNodes_Vel-1));
    U2(Eigen::indexing::all,1) = uNodes(Eigen::seq(nNodes_Vel,2*nNodes_Vel-1));

    Eigen::Vector2d u;
    u = Nu*U2;
    u(0) = abs(u(0));
    u(1) = abs(u(1));

    double magnU = std::sqrt(u.transpose()*u);

    u(0) = magnU; //temporary to force uniform diffusivity
    u(1) = magnU;

    double Pe_Target = 1.0;
    rho = std::max(rho,910.0);
    Eigen::Vector2d diff_art;
    if (magnU == 0){
        diff_art.setZero();
    } else {
        diff_art = u/magnU  *  rho*magnU*scale/2.0/Pe_Target;
    }

    if (Diff_Stab(0,0)<diff_art(0)){
        Diff_Stab(0,0)=diff_art(0);
    }
    if (Diff_Stab(1,1)<diff_art(1)){
        Diff_Stab(1,1)=diff_art(1);
    }

    return Diff_Stab;
}