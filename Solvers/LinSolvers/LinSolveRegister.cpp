#include <iostream>
#include <vector>

#include "LinSolveRegister.h"
#include "MUMPSLinSolver.h"
#include "PARDISOLinSolver.h"
#include "GMRESLinSolver.h"

/// @brief Generates a linear solver based on the requested type
/// @param SolverType Name of the linear solver to be generated
/// @return pointer to linear solver (Remember to delete when not needed)
BaseLinSolver* CreateLinSolver(std::string SolverType){
    std::vector<std::string> PotentialNames = {"BaseLinSolver",
                                               "MUMPS",
                                               "Pardiso",
                                               "GMRES"};
    BaseLinSolver* solver;
    if (SolverType == PotentialNames[0]){
        solver = new BaseLinSolver();
    } else if (SolverType == PotentialNames[1]) {
        solver = new MUMPSLinSolver();
    } else if (SolverType == PotentialNames[2]) {
        solver = new PARDISOLinSolver();
    } else if (SolverType == PotentialNames[3]){
        solver = new GMRESLinSolver();
    } else {
        std::string validTypes = "";
        for (size_t i = 0; i < PotentialNames.size(); i++){
            validTypes.append(PotentialNames[i]);
            validTypes.append(", ");
        }
        throw std::invalid_argument("Linear solver type not defined, valid options are: " + validTypes);
        return 0;
    } 

    return solver;
}