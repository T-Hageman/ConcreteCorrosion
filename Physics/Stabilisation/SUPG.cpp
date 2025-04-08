#include "SUPG.h"

SUPG::SUPG(inputData &inputs){
}

SUPG::~SUPG(){
}

/// @brief Provides streamline-upwind-petrov-galerkin stabilised test functions
/// @param N Non-stabilised test function
/// @param G Non-stabilised test function gradients
/// @param G2 Non-stabilised test function second gradients
/// @param Nu Velocity shape function
/// @param Gu Velocity shape function gradient
/// @param uNodes Nodal velocity vector (u_1; u_2 .. u_n; v_1, v_2 ... v_n)
/// @param scale Typical element lengthscale
/// @param diff Diffusivity equivalent
/// @param N_SUPG output: Stabilised shape function
/// @param dNSUPG_dU output: Stabilised shape function gradient
void SUPG::GetStabilised(   Eigen::RowVectorXd& N, Eigen::MatrixXd& G, Eigen::MatrixXd &G2,
                            Eigen::RowVectorXd &Nu, Eigen::MatrixXd &Gu, Eigen::VectorXd& uNodes, 
                            double scale, double diff, Eigen::RowVectorXd& N_SUPG, Eigen::MatrixXd& dNSUPG_dU){
    // following: Hughes, T. J. R., Scovazzi, G., & Tezduyar, T. E. (2010). Stabilized Methods for Compressible Flows. Journal of Scientific Computing, 43(3), 343–368. https://doi.org/10.1007/s10915-008-9233-5
    //            Eq. 14-17

    double diffLimited = diff;
    if (diffLimited<1.0e-12){
        diffLimited = 1.0e-12;
    }

    size_t nNodes_Vel = uNodes.size()/2;
    Eigen::MatrixXd Nvxvy(2,2*nNodes_Vel);
    Nvxvy.setZero(); Nvxvy(0,Eigen::seq(0,nNodes_Vel-1)) = Nu; Nvxvy(1,Eigen::seq(nNodes_Vel,2*nNodes_Vel-1)) = Nu;

    Eigen::Vector2d u;
    u = Nvxvy*uNodes;
    double uNorm = std::sqrt(u.transpose()*u);
    //bool uZero = false;
    if (uNorm<1.0e-12){
        uNorm = 1.0e-12;
        //uZero = true;
    }

    //Eigen::Matrix2d Grad_u, grad_uNormu;
    //Grad_u = Gu*U2;
    //grad_uNormu = (Eigen::Matrix2d::Identity()/uNorm-u*u.transpose()*std::pow(uNorm, -1.5))*Grad_u;

    //Eigen::MatrixXd G2Mat(2,N.size());
    //for (int i = 0; i < N.size(); i++){
    //    G2Mat(0,i) = u(0)*G2(0,i) + u(1)*G2(2,i);
    //    G2Mat(1,i) = u(0)*G2(2,i) + u(1)*G2(1,i);
    //}
    
    //double Peclet = uNorm*scale/(2*diffLimited); //dimensionless, scale=m (characteristic element size)
    double preFac = 0.5 * scale;// * (1.0/std::tanh(Peclet)-1.0/Peclet) * (1.0-uZero);  // units [m]

    Eigen::Matrix2d dU = Eigen::Matrix2d::Identity()/uNorm-u*u.transpose()/std::pow(uNorm,3);

    N_SUPG = preFac * u.transpose()/uNorm*G;
    dNSUPG_dU = preFac * Nvxvy.transpose()*dU.transpose()*G;
}

/// @brief Generates a representative length scale based on integration point weights
/// @param w Weights of integration points
/// @return scale as used in supg alorithm
double SUPG::getScale(std::vector<double> &w){
    double scale;
    double wSum = std::accumulate(w.begin(), w.end(), 0.0);
    scale = std::sqrt(wSum);
    return scale;
}
