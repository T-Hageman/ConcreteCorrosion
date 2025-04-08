#ifndef SUPG_H
#define SUPG_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../../InputsOutputs/inputData.h"
#include "../../InputsOutputs/SaveData.h"

class SUPG {
    public:
        SUPG(inputData& inputs);
        ~SUPG();

        void GetStabilised(Eigen::RowVectorXd& N, Eigen::MatrixXd& G, Eigen::MatrixXd &G2,
                           Eigen::RowVectorXd &Nu, Eigen::MatrixXd &Gu, Eigen::VectorXd& uNodes, 
                           double scale, double diff, Eigen::RowVectorXd& N_SUPG, Eigen::MatrixXd& G_SUPG);

        double getScale(std::vector<double>& w);
    private:

};

#endif