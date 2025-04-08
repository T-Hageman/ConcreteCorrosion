#ifndef WEAKAREAMODEL_H
#define WEAKAREAMODEL_H

#include <vector>
#include <iostream>
#include <Eigen/Dense>

#include "../BaseModel/BaseModel.h"

class WeakAreaModel: public BaseModel {
    public:
        std::string Name;

        WeakAreaModel(Physics& My_Physics, std::string MyName);
        ~WeakAreaModel();

        void load(inputData& inputs, SaveDataFile& data);
        void save(SaveDataFile& data);

        void init(inputData& inputs);
        void Setup(inputData& inputs);
        void Assemble(Mat& K, Vec& f, Constrainer* cons, size_t step);

        size_t hasTimeData(std::vector<std::string>& DataNames);
        void GetTimeData(std::vector<double>& DataNames);

    protected:

    private:
        size_t ElementGroupIndex;
        std::vector<std::string> dofNames;
        std::vector<size_t> dofSteps;
        std::vector<size_t> dofTypes;

        double xmin, xmax;

        std::vector<double> Vals_0;
        std::vector<double> Vals_t;
        std::vector<double> ExtForce;
        double UpperTime, LowerTime;
        double kdummy;
        bool Lumped;
};

void Register_WeakAreaModel();
BaseModel* New_WeakAreaModel(Physics& My_Physics, std::string MyNameIn);

#endif

