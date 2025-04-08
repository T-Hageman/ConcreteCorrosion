#ifndef B1_4_H
#define B1_4_H

#include <vector>
#include <Eigen/Dense>

#include "../BaseElemType.h"
#include "../IPScheme.h"

/// @brief Quartic Bezier line element
class B1_4: public BaseElemType {
    public:
        B1_4(int ipcount1D);
        ~B1_4();

        void getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData);
        //void getShapeGrads2(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData);
        void getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData);
        void getNormals(std::vector<Eigen::VectorXd> &normals, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData);

    private:
        void SetUpBaseShapes();

};

#endif