#ifndef T2_6B_H
#define T2_6B_H

#include <vector>
#include <Eigen/Dense>

#include "../BaseElemType.h"
#include "../IPScheme.h"

/// @brief Quartic Bezier line element
class T2_6B: public BaseElemType {
    public:
        T2_6B(int ipcount1D);
        ~T2_6B();

        void getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes);

    private:
        void SetUpBaseShapes();

};

#endif