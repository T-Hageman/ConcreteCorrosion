#ifndef L4_H
#define L4_H

#include <vector>
#include <Eigen/Dense>

#include "../BaseElemType.h"
#include "../IPScheme.h"

/// @brief 4 Node square Lagrangian plane element
class L4: public BaseElemType {
    public:
        L4(int ipcount1D);
        ~L4();

        void getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes);

    private:
        void SetUpBaseShapes();

};

#endif