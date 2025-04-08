#ifndef PlotData_XY_H
#define PlotData_XY_H
#ifdef LOCALENVIRONMENT
#include "Visualisation.h"

class PlotData_XY: public Vis_Plot{

    public:
        PlotData_XY(std::string name, inputData& inputs, Physics* phys);
        ~PlotData_XY();
        void Plot();



    private:
        std::string xTypeString;
        std::vector<std::string> yTypeString;
        std::vector<size_t> typeIndices;
        vtkNew<vtkChartXY> chart;

        vtkNew<vtkTable> PlotXY_table;
        std::vector<vtkPlot*> lines;
        std::string xlim, ylim;

};

#endif
#endif