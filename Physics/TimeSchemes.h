#ifndef TIMESCHEMES_H
#define TIMESCHEMES_H

#include <vector>
#include <iostream>

#include "../InputsOutputs/inputData.h"
#include "../utility/utility.h"

/// @brief Time discretisation schemes implmentation
class TimeScheme {
    public:
        TimeScheme(inputData& inputs);
        ~TimeScheme();
        void UpdateVelAcc(GhostVector& State, GhostVector& dState, GhostVector& ddState, GhostVector& StateOld, GhostVector& dStateOld, GhostVector& ddStateOld);
        void Set_dt(double dt_in);

        double dt;      //time increment
        double du_dt;   //derivate of state change to state (e.g. dv/dx)
        double ddu_dt;  //second derivate of state change to state (e.g. da/dx)
    private:
        std::vector<double> TimeConst;  //constants used within time discretisation
        std::string DiscType;           //type of time discretisation
};

#endif