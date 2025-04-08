#include "Visualisation.h"
#ifdef LOCALENVIRONMENT

#include <iostream>
#include <set>

#include "../physics.h"
#include "Plot_Surf2D.h"
#include "Plot_Surf3D.h"
#include "PlotData_XY.h"

#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCylinderSource.h>
#include <vtkNew.h>
#include <vtkPolyDataMapper.h>
#include <vtkCellType.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkChartXYZ.h>
#include <vtkPlotSurface.h>
#include <vtkTable.h>
#include <vtkChartXYZ.h>
#include <vtkChartXY.h>
#include <vtkContextMouseEvent.h>
#include <vtkContextScene.h>
#include <vtkContextView.h>
#include <vtkFloatArray.h>
#include <vtkNamedColors.h>
#include <vtkLookupTable.h>
#include <vtkNew.h>
#include <vtkPen.h>
#include <vtkPlot.h>
#include <vtkPlotSurface.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkTable.h>
#include <vtkVersion.h>
#include <vtkCubeAxesActor.h>
#include <vtkColorSeries.h>
#include <vtkVariantArray.h>
#include <vtkAxis.h>
#include <vtkChartLegend.h>
#include <vtkTextActor.h>
#include <vtkTextProperty.h>
#include <vtkDataArray.h>


Visualisation::Visualisation(inputData& inputs, Physics* phys){
    if (inputs.HasKey({"Visualisation"})){
        std::vector<std::string> PlotNames; inputs.GetRequired(PlotNames, {"Visualisation","Plots"});
        nPlots = PlotNames.size();
        All_Plots.resize(nPlots); 
        for (size_t i = 0; i < nPlots; i++){
            std::string PlotTypeName;
            PlotTypes PlotType;
            Logs.PrintSingle("Plotting: "+PlotNames[i]+"\n",3);
            inputs.GetRequired(PlotTypeName, {"Visualisation",PlotNames[i],"Type"});
            if (PlotTypeMap.count(PlotTypeName)){
                PlotType = PlotTypeMap[PlotTypeName];
            } else {
                throw std::invalid_argument("Visualisation does not recognise plot type "+PlotTypeName);
            }

            switch (PlotType) {
            case PlotTypes::Mesh_2D:
                All_Plots[i] = new Vis_Plot(PlotNames[i], inputs, phys);
                break;
            case PlotTypes::Surf_2D:
                All_Plots[i] = new Plot_Surf2D(PlotNames[i], inputs, phys);
                break;
            case PlotTypes::Surf_3D:
                All_Plots[i] = new Plot_Surf3D(PlotNames[i], inputs, phys);
                break;
            case PlotTypes::Plot_XY:
                All_Plots[i] = new PlotData_XY(PlotNames[i], inputs, phys);
                break;
            }
        }
    } else {
        nPlots = 0;
    }

    Plot("once", 0);
}

Visualisation::~Visualisation(){

}

void Visualisation::Plot(std::string step, size_t stp){
    for (size_t i = 0; i < nPlots; i++){
        if (step == "once"){
            All_Plots[i]->Plot();
        }
        if (step == All_Plots[i]->PlotFrequency){
            if (All_Plots[i]->PlotFrequency == "step" || (All_Plots[i]->PlotFrequency == "it" && stp == All_Plots[i]->DofStep_Data)){
                MPI_Barrier(PETSC_COMM_WORLD);
                Logs.PrintSingle("Updating Figure:"+All_Plots[i]->Name+"\n",3);
                MPI_Barrier(PETSC_COMM_WORLD);
                All_Plots[i]->Plot();
                MPI_Barrier(PETSC_COMM_WORLD);
                Logs.PrintSingle("Finished updating Figure:"+All_Plots[i]->Name+"\n",3);
                MPI_Barrier(PETSC_COMM_WORLD);
            }
        }
    }
}


Vis_Plot::Vis_Plot(std::string name, inputData& inputs, Physics* phys){
    PetscCallThrow(MPI_Comm_size(PETSC_COMM_WORLD,&size));
    PetscCallThrow(MPI_Comm_rank(PETSC_COMM_WORLD,&rank));
    DoOnce = true;

    Name = name;
    physics = phys;

    inputs.GetRequired(PlotFrequency, {"Visualisation",Name,"Frequency"});

    renderWindow = view->GetRenderWindow();
    if (rank==0){
        renderWindow->SetSize(600, 600);
        renderWindow->Render();
    }
}

Vis_Plot::~Vis_Plot(){

}

void Vis_Plot::Plot(){

}

void Vis_Plot::AddColourBar(){
    std::stringstream colorbarTitle;
    if (LogScale){
        colorbarTitle << "log10(" << DataName << ")";
    } else {
        colorbarTitle << DataName;
    }
    ColourBar->SetTitle(colorbarTitle.str().c_str());
    ColourBar->GetTitleTextProperty()->SetColor(
        colors->GetColor3d("black").GetData());
    ColourBar->GetLabelTextProperty()->SetColor(
        colors->GetColor3d("black").GetData());
    ColourBar->GetAnnotationTextProperty()->SetColor(
        colors->GetColor3d("black").GetData());
    //ColourBar->UnconstrainedFontSizeOn();
    ColourBar->SetNumberOfLabels(5);
    ColourBar->SetOrientation(0);
    if (LogScale){
        ColourBar->SetLabelFormat("%.2f");
    } else {
        ColourBar->SetLabelFormat("%4.3e");
    }
    ColourBar->GetTitleTextProperty()->SetBold(false);
    ColourBar->GetTitleTextProperty()->SetItalic(false);
    ColourBar->GetLabelTextProperty()->SetBold(false);
    ColourBar->GetLabelTextProperty()->SetItalic(false);

    ColourBar->SetWidth(500);
    ColourBar->SetMaximumWidthInPixels(500);
    ColourBar->SetMaximumHeightInPixels(40);
    ColourBar->SetPosition(0.1, 0.025);

    renderer->AddActor2D(ColourBar);
}

void Vis_Plot::AddGrid(double GridRanges[6]){
    vtkNew<vtkCubeAxesActor> cubeAxesActor;
    vtkColor3d axis1Color = colors->GetColor3d("Salmon");
    vtkColor3d axis2Color = colors->GetColor3d("PaleGreen");
    vtkColor3d axis3Color = colors->GetColor3d("LightSkyBlue");

    cubeAxesActor->SetBounds(GridRanges);
    cubeAxesActor->SetCamera(renderer->GetActiveCamera());
    cubeAxesActor->GetTitleTextProperty(0)->SetColor(axis1Color.GetData());
    cubeAxesActor->GetTitleTextProperty(0)->SetFontSize(48);
    cubeAxesActor->GetLabelTextProperty(0)->SetColor(axis1Color.GetData());
    cubeAxesActor->SetXTitle("x [m]");
    cubeAxesActor->SetYTitle("y [m]");
    cubeAxesActor->SetZTitle("z [m]");

    cubeAxesActor->GetTitleTextProperty(1)->SetColor(axis2Color.GetData());
    cubeAxesActor->GetLabelTextProperty(1)->SetColor(axis2Color.GetData());

    cubeAxesActor->GetTitleTextProperty(2)->SetColor(axis3Color.GetData());
    cubeAxesActor->GetLabelTextProperty(2)->SetColor(axis3Color.GetData());

    cubeAxesActor->DrawXGridlinesOff();
    cubeAxesActor->DrawYGridlinesOff();
    cubeAxesActor->DrawZGridlinesOff();
    cubeAxesActor->SetGridLineLocation(cubeAxesActor->VTK_GRID_LINES_FURTHEST);

    cubeAxesActor->XAxisMinorTickVisibilityOff();
    cubeAxesActor->YAxisMinorTickVisibilityOff();
    cubeAxesActor->ZAxisMinorTickVisibilityOff();

    cubeAxesActor->SetFlyModeToStaticEdges();

    renderer->AddActor(cubeAxesActor);
}






#endif