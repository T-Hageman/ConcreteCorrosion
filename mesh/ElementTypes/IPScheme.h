#ifndef IPSCHEME_H
#define IPSCHEME_H

#include <vector>
#include <Eigen/Dense>

void getIPScheme1D(std::vector<double> *wout, std::vector<double> *coordsIPout, int ipCount);
int getIPScheme(std::vector<double> *wout, std::vector<Eigen::VectorXd> *coordsIPout, int ipCount, int dim);
int getIPSchemeTriangle(std::vector<double> *wout, std::vector<Eigen::VectorXd> *coordsIPout, int ipCount, int dim);

#endif