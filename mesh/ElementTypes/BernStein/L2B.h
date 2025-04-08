#ifndef L2B_H
#define L2B_H

#include <vector>
#include <Eigen/Dense>

#include "../BaseElemType.h"
#include "../IPScheme.h"

/// @brief 2 Node Bezier line element
class L2B: public BaseElemType {
    public:
        L2B(int ipcount1D);
        ~L2B();

        void getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes);

    private:
        void SetUpBaseShapes();

};

#endif