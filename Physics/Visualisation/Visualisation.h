#ifdef LOCALENVIRONMENT
#ifndef VISUALISATION_H
#define VISUALISATION_H

#include "../../InputsOutputs/inputData.h"

#include <string>
#include <vector>
#include <vtkRenderer.h>
#include <vtkNamedColors.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkPointData.h>
#include <vtkScalarBarActor.h>
#include <vtkDataSetMapper.h>
#include <vtkTable.h>
#include <vtkChartXY.h>

class Physics;

class Vis_Plot {
    public:
        Vis_Plot(std::string name, inputData& inputs, Physics* phys);
        virtual ~Vis_Plot();
        virtual void Plot();

        std::string PlotFrequency;
        size_t DofStep_Data;

        std::string Name;
    protected:
        PetscMPIInt    size;
        PetscMPIInt    rank;

        Physics* physics;
        bool DoOnce;

        
        //graphics objects
        vtkNew<vtkRenderer> renderer;
        vtkNew<vtkContextView> view;
        vtkRenderWindow* renderWindow;
        vtkNew<vtkNamedColors> colors;
        vtkNew<vtkTextActor> textActor;

        void AddGrid(double GridRanges[6]);
        void AddColourBar();

        size_t ElemGroupIndex_Data;

        std::string DataName;

        bool LogScale;

        //display mesh data
        vtkNew<vtkPoints> points;
        vtkNew<vtkUnstructuredGrid> ugrid;
        vtkNew<vtkDataSetMapper> ugridMapper;
        vtkNew<vtkActor> ugridActor;
        vtkNew<vtkScalarBarActor> ColourBar;

        private:
        
};

/// @brief Plots figures while running
class Visualisation {
    public:
        Visualisation(inputData& inputs, Physics* phys);
        ~Visualisation();

        void Plot(std::string step, size_t stp);

        enum PlotTypes{Mesh_2D, Surf_2D, Surf_3D, Plot_XY, BoundaryFill, Contour3D, IPSurf_3D};
        std::unordered_map<std::string, PlotTypes> PlotTypeMap = {
            {"Mesh_2D", Mesh_2D},
            {"Surf_2D", Surf_2D},
            {"Surf_3D", Surf_3D},
            {"Plot_XY", Plot_XY}
        };
    private:
        size_t nPlots;
        std::vector<Vis_Plot*> All_Plots;
};


#endif
#endif