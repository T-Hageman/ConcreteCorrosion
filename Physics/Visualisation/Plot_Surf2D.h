#ifndef PLOT_SURF2D_H
#define PLOT_Surf2D_H
#ifdef LOCALENVIRONMENT
#include "Visualisation.h"

class Plot_Surf2D: public Vis_Plot{

    public:
        Plot_Surf2D(std::string name, inputData& inputs, Physics* phys);
        ~Plot_Surf2D();
        void Plot();



    private:
        std::vector<std::string> DofNames_u = {"ux","uy","uz"};
        size_t dofStep_u;
        std::vector<size_t> dofTypes_u;

        size_t ElemGroupIndex_u;

        bool Deformations;
        double DeformationScale;

        bool ShowMesh;

        void DataToUGrid(std::vector<std::vector<double>> &XCoords, std::vector<std::vector<double>> &YCoords, std::vector<std::vector<double>> &Data, 
                std::vector<std::vector<double>> &UX, std::vector<std::vector<double>> &UY, bool HasColour, vtkNew<vtkLookupTable>& colorLookupTable, bool UpdateCOnly);

};

#endif
#endif