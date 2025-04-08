#ifndef T3_10_H
#define T3_10_H

#include <vector>
#include <Eigen/Dense>

#include "../BaseElemType.h"
#include "../IPScheme.h"

/// @brief Cubic Bezier plane element
class T3_10: public BaseElemType {
    public:
        T3_10(int ipcount1D);
        ~T3_10();

        void getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes);

    private:
        void SetUpBaseShapes();

};

#endif