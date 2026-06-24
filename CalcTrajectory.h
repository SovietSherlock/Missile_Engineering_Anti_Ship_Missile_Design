#ifndef CALCTRAJECTORY_H
#define CALCTRAJECTORY_H

#include <vector>
#include <string>
#include "DataStructs.h"
#include "LibConstFunc.h"

std::string CalcTrajectory(
    const CalcParameters      &calc,
    const InitConditions      &inits,
    const LimitConditions     &limits,
    const C_Data              &Cdata,
    const N_Data              &Ndata,
    const R_Data              &Rdata,
    const GuidanceMethod_Data &MethodData,
    const RocketParams        &rocket,
    double mu_0,
    double eta_0,
    double K_m,
    double K_P,
    double K_g,
    std::vector<TrajectoryParameters> &s);

#endif // CALCTRAJECTORY_H
