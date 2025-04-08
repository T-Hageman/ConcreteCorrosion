#ifndef PLOT_SURF3D_H
#define PLOT_SURF3D_H
#ifdef LOCALENVIRONMENT
#include "Visualisation.h"

class Plot_Surf3D: public Vis_Plot{

    public:
        Plot_Surf3D(std::string name, inputData& inputs, Physics* phys);
        ~Plot_Surf3D();
        void Plot();

    private:
        std::vector<std::string> DofNames_u = {"ux","uy","uz"};

        size_t ElemGroupIndex_u;

        void DataToUGrid(std::vector<std::vector<double>> &XCoords, std::vector<std::vector<double>> &YCoords, std::vector<std::vector<double>> &ZCoords, 
                          std::vector<std::vector<double>> &Data, std::vector<std::vector<double>> &UX, std::vector<std::vector<double>> &UY, std::vector<std::vector<double>> &UZ,
                          vtkNew<vtkLookupTable>& colorLookupTable, bool UpdateCOnly);

        bool Deformations;
        double DeformationScale;
};

#endif
#endif