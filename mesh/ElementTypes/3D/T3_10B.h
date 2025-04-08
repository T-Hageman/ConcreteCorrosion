#ifndef T3_10B_H
#define T3_10B_H

#include <vector>
#include <Eigen/Dense>

#include "../BaseElemType.h"
#include "../IPScheme.h"

/// @brief Cubic Bezier plane element
class T3_10B: public BaseElemType {
    public:
        T3_10B(int ipcount1D);
        ~T3_10B();

        void getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes);

    private:
        void SetUpBaseShapes();

};

#endif