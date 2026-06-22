#ifndef CALCTRAJECTORY_H
#define CALCTRAJECTORY_H

#include <vector>
#include <string>
#include "DataStructs.h"

enum class TrajectoryMode {
    EWT,      // Энергетически выгодная траектория
    GORKA     // С манёвром «горка» на конечном участке
};

std::string CalcTrajectory_Euler(
    const CalcParameters& calcParams,
    const InitConditions& inits,
    const LimitConditions& limits,
    const C_Data& cdata,
    const N_Data& ndata,
    const R_Data& rdata,
    const GuidanceMethod_Data& methodData,
    TrajectoryMode mode,
    std::vector<TrajectoryParameters>& result);

#endif // CALCTRAJECTORY_H