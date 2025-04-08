#include "PlotData_XY.h"
#ifdef LOCALENVIRONMENT

#include <iostream>
#include <set>

#include "../physics.h"

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



PlotData_XY::PlotData_XY(std::string name, inputData &inputs, Physics *phys): Vis_Plot(name, inputs, phys) {
    if (rank == 0){
        bool HasDataKey;

        inputs.GetRequired(xTypeString, {"Visualisation",Name,"x"});
        inputs.GetRequired(yTypeString, {"Visualisation",Name,"y"});
        typeIndices.resize(yTypeString.size()+1);
        
        vtkNew<vtkFloatArray> arrX;
        arrX->SetName(xTypeString.c_str());
        PlotXY_table->AddColumn(arrX);
        typeIndices[0] = Find(physics->TimeDataTypes, xTypeString, HasDataKey);
        if (HasDataKey==false){
            throw std::invalid_argument("Dataplot has no "+xTypeString+"\n");
        }

        std::vector<vtkNew<vtkFloatArray>> arrsY(yTypeString.size());
        for (size_t i = 0; i < yTypeString.size(); i++){
            arrsY[i]->SetName(yTypeString[i].c_str());
            PlotXY_table->AddColumn(arrsY[i]);
            typeIndices[i+1] = Find(physics->TimeDataTypes, yTypeString[i], HasDataKey);
            if (HasDataKey==false){
                throw std::invalid_argument("Dataplot has no "+xTypeString+"\n");
            }
        }  
        PlotXY_table->InsertNextBlankRow(1e-12);

        xlim = ""; ylim = "";
        inputs.GetOptional(ylim, {"Visualisation",Name,"ylim"});
        inputs.GetOptional(xlim, {"Visualisation",Name,"xlim"});
    }
}

PlotData_XY::~PlotData_XY(){

}

void PlotData_XY::Plot(){
    if (DoOnce){
        if (rank == 0){
            view->GetScene()->AddItem(chart);

            vtkNew<vtkColorSeries> colorSeries;
            colorSeries->SetColorScheme(vtkColorSeries::SPECTRUM);

            if (ylim=="log"){
                chart->GetAxis(vtkAxis::LEFT)->LogScaleOn();
                chart->GetAxis(vtkAxis::LEFT)->SetLogScale(true);
                chart->GetAxis(vtkAxis::LEFT)->SetMinimum(1e-3);
                chart->GetAxis(vtkAxis::LEFT)->SetMinimumLimit(1e-9);
                chart->GetAxis(vtkAxis::LEFT)->Update();
            }

            lines.resize(typeIndices.size()-1);
            for (size_t i = 0; i < typeIndices.size()-1; i++){
                lines[i] = chart->AddPlot(vtkChart::LINE);
                auto c = colorSeries->GetColorRepeating(i);
                lines[i]->SetColorF(c[0], c[1], c[2]);
                lines[i]->SetWidth(3.0);
                lines[i]->SetInputData(PlotXY_table, 0, i+1);
            }
            chart->GetAxis(vtkAxis::BOTTOM)->SetTitle(xTypeString.c_str());
            chart->SetShowLegend(true);

            vtkChartLegend* Leg=chart->GetLegend();
            Leg->SetHorizontalAlignment(Leg->LEFT);
        }
        DoOnce = false;
    }
    if (rank == 0){
        bool x_above_zero = false;
        size_t rowIdx = PlotXY_table->InsertNextBlankRow();
        for (size_t i = 0; i < typeIndices.size(); i++){
            size_t lastIndex = physics->TimeData[typeIndices[i]].size()-1;
            double pdata = physics->TimeData[typeIndices[i]][lastIndex];
            PlotXY_table->SetValue(rowIdx, i, pdata);
            //std::cout << lastIndex+1 << ", " << i << ",  " << pdata << "\n";
            if (i==0 && pdata>0.0){
                x_above_zero = true;
            }
        }
        PlotXY_table->Modified();
        
        for (size_t i = 0; i < typeIndices.size()-1; i++){
            lines[i]->Update();
        }

        chart->Update();
        chart->RecalculateBounds();

	    double XMin = chart->GetAxis(vtkAxis::BOTTOM)->GetMinimum();
        double XMax = chart->GetAxis(vtkAxis::BOTTOM)->GetMaximum();
        if (xlim=="0+" && x_above_zero){
            chart->GetAxis(vtkAxis::BOTTOM)->SetMinimum(0.0);
            chart->GetAxis(vtkAxis::BOTTOM)->SetMinimumLimit(0.0);
            chart->GetAxis(vtkAxis::BOTTOM)->Update();
        }

        //Logatithmic scale
        if (ylim=="log"){
            //chart->GetAxis(vtkAxis::LEFT)->LogScaleOn();
            //chart->GetAxis(vtkAxis::LEFT)->SetLogScale(true);
            //chart->GetAxis(vtkAxis::LEFT)->SetMinimum(1e-3);
            //chart->GetAxis(vtkAxis::LEFT)->SetMinimumLimit(1e-3);
            //chart->GetAxis(vtkAxis::LEFT)->Update();
        }

        renderWindow->Render();
    }
}
#endif
