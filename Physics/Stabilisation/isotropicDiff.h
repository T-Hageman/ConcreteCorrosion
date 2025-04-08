#ifndef ISOTROPICDIFF_H
#define ISOTROPICDIFF_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../../InputsOutputs/inputData.h"
#include "../../InputsOutputs/SaveData.h"

double GetIsotropicDiffusivity(Eigen::RowVectorXd &Nu,Eigen::VectorXd& uNodes, std::vector<double>& w, double diff);
Eigen::Matrix4d GetIsotropicDiffusivity(Eigen::RowVectorXd &Nu,Eigen::VectorXd& uNodes, std::vector<double>& w, Eigen::Matrix4d& diff);
Eigen::Matrix4d GetIsotropicDiffusivity(Eigen::RowVectorXd &Nu,Eigen::VectorXd& uNodes, std::vector<double>& w, Eigen::Matrix4d& diff, double rho);

#endif