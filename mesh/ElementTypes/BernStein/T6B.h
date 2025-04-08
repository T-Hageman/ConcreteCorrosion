#ifndef T6B_H
#define T6B_H

#include <vector>
#include <Eigen/Dense>

#include "../BaseElemType.h"
#include "../IPScheme.h"

/// @brief 6 Node triangular Bernstein plane element
class T6B: public BaseElemType {
    public:
        T6B(int ipcount1D);
        ~T6B();

        void getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes);

    private:
        void SetUpBaseShapes();

};

#endif