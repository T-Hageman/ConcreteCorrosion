#ifndef B1_3_H
#define B1_3_H

#include <vector>
#include <Eigen/Dense>

#include "../BaseElemType.h"
#include "../IPScheme.h"

/// @brief Cubic Bezier line element
class B1_3: public BaseElemType {
    public:
        B1_3(int ipcount1D);
        ~B1_3();

        void getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData);
        //void getShapeGrads2(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData);
        void getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData);
        void getNormals(std::vector<Eigen::VectorXd> &normals, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData);
        
    private:
        void SetUpBaseShapes();

};

#endif