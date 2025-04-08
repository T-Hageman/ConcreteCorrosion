#ifndef BASEELEMTYPE_H
#define BASEELEMTYPE_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "IPScheme.h"

void BezierBasis(int p, double x_ip, Eigen::VectorXd& N, Eigen::VectorXd& G, Eigen::VectorXd& G2);
int fact(int n);

/// @brief Parametric element (specific implementation performed inside derived classes)
class BaseElemType {
    public:
        std::string Name;

        BaseElemType();
        virtual ~BaseElemType();

        void init(int ipcount, int dim);
        virtual void getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes);
        virtual void getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData);
        virtual void getShapeGrads(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData);

        virtual void getShapeGrads2(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes);
        virtual void getShapeGrads2(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData);
        virtual void getShapeGrads2(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<Eigen::MatrixXd> &G2out, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData);

        virtual void getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes);
        virtual void getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData);
        virtual void getCoordsIP(Eigen::MatrixXd &coordsIP, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData);

        virtual void getNormals(std::vector<Eigen::VectorXd> &normals, Eigen::MatrixXd &coordsNodes);
        virtual void getNormals(std::vector<Eigen::VectorXd> &normals, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData);
        virtual void getNormals(std::vector<Eigen::VectorXd> &normals, Eigen::MatrixXd &coordsNodes, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData);

        virtual void getExportShape(std::vector<Eigen::RowVectorXd> &Nout);
        virtual void getExportShape(std::vector<Eigen::RowVectorXd> &Nout, Eigen::MatrixXd &ElemData);
        virtual void getExportShape(std::vector<Eigen::RowVectorXd> &Nout, Eigen::MatrixXd &ElemData, Eigen::VectorXd &NodeData);

        uint NodeCount; //Number of nodes associated with the element
        uint dim;   //spatial dimension this element exists in
        uint ipcount;   //Number of integration points

        bool requiresData;  //Flag indication if elemental data is required for this element
        bool requiresNodeData;  //Flag indication if nodal data is required for this element

        std::vector<double> w;  //integration point weights for parametric element
        std::vector<Eigen::VectorXd> coordsIP;  //parametric coordinates of integration points

        std::vector<Eigen::RowVectorXd> Nbase, NExport; //Shape functions, N[ip] to obtain the rowvector of shape functions within integration point
        std::vector<Eigen::MatrixXd> Gbase; //shape function gradients in parametric space, G[ip] to obtain matrix (x/y x N)
        std::vector<Eigen::MatrixXd> G2base;    //Second gradients of shape functions in parametric space G2[ip] to obtain matrix (xx/yy/xy x N)

        std::vector<size_t> PlottingOrder; //Used for visualisation, order points are provided in for polygon
    protected:
        void getShapeGrads_2D(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes);
        void getShapeGrads_2D_Line(std::vector<Eigen::RowVectorXd> &Nout, std::vector<Eigen::MatrixXd> &Gout, std::vector<double> &wout, Eigen::MatrixXd &coordsNodes);

    private:

};

#endif

