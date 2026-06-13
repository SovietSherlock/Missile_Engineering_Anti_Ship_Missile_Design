#ifndef CALCTRAJECTORY_H
#define CALCTRAJECTORY_H

#include <vector>
#include <string>

#include "DataStructs.h"

std::string CalcTrajectory_Euler(const CalcParameters& CalcParameters,
    const InitConditions& inits,
    const LimitConditions& limits,
    const C_Data& Cdata,
    const N_Data& Ndata,
    const R_Data& Rdata,
    const GuidanceMethod_Data& MethodData,
    std::vector<TrajectoryParameters>& result);


#endif // CALCTRAJECTORY_H
