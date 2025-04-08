#ifndef AREACONSTRAINTSMODEL_H
#define AREACONSTRAINTSMODEL_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"

class AreaConstraintsModel: public BaseModel {
    public:
        std::string Name;

        AreaConstraintsModel(Physics& My_Physics, std::string MyName);
        ~AreaConstraintsModel();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);
    protected:

    private:
        size_t ElementGroupIndex;
        std::vector<std::string> dofNames;
        std::vector<size_t> dofSteps;
        std::vector<size_t> dofTypes;
        std::vector<double> conVals;

        std::vector<double> xLims, yLims, zLims;
        double kDummy;
};

void Register_AreaConstraintsModel();
BaseModel* New_AreaConstraintsModel(Physics& My_Physics, std::string MyNameIn);

#endif

