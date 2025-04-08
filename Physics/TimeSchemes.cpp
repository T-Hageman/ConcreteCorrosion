#include "TimeSchemes.h"

/// @brief Initialize the time discretisation scheme based on inputs provided
/// @param inputs Input parameter object
TimeScheme::TimeScheme(inputData& inputs){
    inputs.GetRequired(DiscType, {"properties","TimeDiscretisation","Scheme"});
    if (DiscType == "Newmark"){
        TimeConst.resize(2);
        inputs.GetRequired(TimeConst[0], {"properties","TimeDiscretisation","beta"});
        inputs.GetRequired(TimeConst[1], {"properties","TimeDiscretisation","gamma"});
    } else if (DiscType == "Alpha") {
        TimeConst.resize(1);
        inputs.GetRequired(TimeConst[0], {"properties","TimeDiscretisation","alpha"});
    } else if (DiscType == "Euler") {
        DiscType = "Alpha";
        TimeConst.resize(1);
        TimeConst[0] = 1.0;
    } else if (DiscType == "BDF2") {

    } else {
        throw std::invalid_argument("Time Discretisation type not defined");
    }
}

TimeScheme::~TimeScheme(){

}

/// @brief Updates the velocity and acceleration based on current and old state vectors
/// @param State x^{t+Delta t}
/// @param dState v^{t+Delta t}
/// @param ddState a^{t+Delta t}
/// @param StateOld x^t
/// @param dStateOld v^t
/// @param ddStateOld a^t
void TimeScheme::UpdateVelAcc(GhostVector& State, GhostVector& dState, GhostVector& ddState, GhostVector& StateOld, GhostVector& dStateOld, GhostVector& ddStateOld){
    if (DiscType == "Newmark"){
        // v = gamma/(beta*dt) * (x-xold) - (gamma/beta-1)*vold - (dt*gamma/(2beta)-dt)*aold
        // a = 1/(beta*dt^2) * (x-xold) - 1/(beta*dt)*vold - (1/(2beta)-1)*aold
        VecAXPBYPCZ(ddState.DataVector, 1.0/(dt*dt*TimeConst[0]), -1.0/(dt*dt*TimeConst[0]), 0.0, State.DataVector, StateOld.DataVector);
        VecAXPBYPCZ(ddState.DataVector, -1.0/(TimeConst[0]*dt), -(1.0/(2.0*TimeConst[0])-1.0), 1.0, dStateOld.DataVector, ddStateOld.DataVector);

        VecAXPBYPCZ(dState.DataVector, TimeConst[1]/(dt*TimeConst[0]), -TimeConst[1]/(dt*TimeConst[0]), 0.0, State.DataVector, StateOld.DataVector);
        VecAXPBYPCZ(dState.DataVector, -(TimeConst[1]/TimeConst[0]-1.0), -(dt*TimeConst[1]/(2.0*TimeConst[0])-dt), 1.0, dStateOld.DataVector, ddStateOld.DataVector);
    } else if (DiscType == "Alpha") {
        // v = 1/(alpha*dt)*(x-xold) + (1-1/alpha)*vold     (Euler for alpha = 1)
        // a = 0
        ddState.Zero();
        VecAYPX(dState.DataVector, 0.0, dStateOld.DataVector);
        VecAXPBYPCZ(dState.DataVector, 1.0/(dt*TimeConst[0]), -1.0/(dt*TimeConst[0]), (1.0-1.0/TimeConst[0]), State.DataVector, StateOld.DataVector);
    } else if (DiscType == "BDF2") { //STILL TO CHECK
        // v = 1.5 (x-xold)/dt - vold/2 + dt/4*aold
        // a = (x-xold)/dt^2 - vold/dt

        VecAXPBYPCZ(ddState.DataVector, 1.0/(dt*dt), -1.0/(dt*dt), 0.0, State.DataVector, StateOld.DataVector);
        VecAXPY(ddState.DataVector, -1.0, dStateOld.DataVector);
        
        VecAXPBYPCZ(dState.DataVector, 1.0/dt, -1.0/dt, 0.0, State.DataVector, StateOld.DataVector);
        VecAXPBYPCZ(dState.DataVector, -0.0, 0.0*dt, 1.0, dStateOld.DataVector, ddStateOld.DataVector);
    }
}

/// @brief Sets time increment size and updates derivatives
/// @param dt_in New time increment size to use
void TimeScheme::Set_dt(double dt_in){
    dt = dt_in;
    if (DiscType == "Newmark"){
        du_dt  = TimeConst[1]/(dt*TimeConst[0]);
        ddu_dt = 1.0/(dt*dt*TimeConst[0]);
    } else if (DiscType == "Alpha") {
        du_dt  = 1.0/(dt*TimeConst[0]);
        ddu_dt = 0.0;
    } else if (DiscType == "BDF2") {
        du_dt  = 1.0/dt;
        ddu_dt = 1.0/(dt*dt); 
    }
}