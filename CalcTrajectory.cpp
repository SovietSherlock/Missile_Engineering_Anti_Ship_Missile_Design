#include "CalcTrajectory.h"
#include <cmath>
#include <cfloat>
#include "LibConstFunc.h"
#include "KinematicFunc2D.h"
#include "GuidanceMethods.h"
#include "Atm_GOST4401.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// ============================================================================
// Вспомогательная функция проверки NaN
// ============================================================================
static inline bool isNaN(double x)
{
    return x != x;
}

// ============================================================================
// 1. КИНЕМАТИЧЕСКАЯ МОДЕЛЬ (без физики)
// ============================================================================
std::string CalcTrajectory_Euler(
    const CalcParameters& calcParams,
    const InitConditions& inits,
    const LimitConditions& limits,
    const C_Data& Cdata,
    const N_Data& Ndata,
    const R_Data& Rdata,
    const GuidanceMethod_Data& MethodData,
    std::vector<TrajectoryParameters>& s)
{
    double dt = calcParams.dt;
    unsigned int N_max = static_cast<unsigned int>(ceil(limits.t_r_max / dt));

    s.clear();
    s.reserve(N_max);

    // === Инициализация ===
    TrajectoryParameters p = {};
    p.t = 0.0;
    p.Theta_c = inits.Theta_c0;
    p.v_c = inits.v_c0;
    p.x_g_c = inits.x_g_c0;
    p.y_g_c = inits.y_g_c0;
    p.Theta_n = inits.Theta_n0;
    p.v_n = inits.v_n0;
    p.x_g_n = inits.x_g_pusk;
    p.y_g_n = inits.y_g_pusk;
    p.Theta_r = inits.Theta_r0;
    p.v_r = inits.v_r0;
    p.x_g_r = inits.x_g_pusk;
    p.y_g_r = inits.y_g_pusk;

    s.push_back(p);

    // === Начальные углы ===
    if (isNaN(s[0].Theta_r))
    {
        s[0].epsilon_rc = epsilon(s[0].x_g_r, s[0].y_g_r, s[0].x_g_c, s[0].y_g_c);
        s[0].epsilon_nc = s[0].epsilon_rc;
        s[0].dot_epsilon_rc = dot_epsilon(s[0].x_g_r, s[0].y_g_r, s[0].x_g_c, s[0].y_g_c,
                                          s[0].v_r, s[0].v_c, s[0].epsilon_rc, s[0].Theta_c);
        s[0].dot_epsilon_nc = dot_epsilon(s[0].x_g_n, s[0].y_g_n, s[0].x_g_c, s[0].y_g_c,
                                          s[0].v_n, s[0].v_c, s[0].Theta_n, s[0].Theta_c);

        switch (MethodData.Method)
        {
        case GuidanceMethod::Null:
        case GuidanceMethod::PurePursuit:
        case GuidanceMethod::ProportionalNav_const:
        case GuidanceMethod::ProportionalNav_var:
        case GuidanceMethod::EVT:
            s[0].Theta_r = s[0].epsilon_rc;
            break;

        case GuidanceMethod::DeviatedPursuit:
            s[0].Theta_r = s[0].epsilon_rc + MethodData.phi_upr * sign(s[0].dot_epsilon_rc);
            break;

        case GuidanceMethod::ParallelNav:
            s[0].Theta_r = s[0].epsilon_rc - asin(s[0].v_c / s[0].v_r * sin(s[0].epsilon_rc - s[0].Theta_c));
            break;

        case GuidanceMethod::BeamRiding_LineOfSight:
            s[0].epsilon_l = s[0].epsilon_nc;
            s[0].Theta_r = s[0].epsilon_l + asin(s[0].v_n / s[0].v_r * sin(s[0].Theta_n - s[0].epsilon_l));
            break;

        case GuidanceMethod::BeamRiding_ConstAlignment:
            s[0].Delta_r = r(s[0].x_g_n, s[0].y_g_n, s[0].x_g_c, s[0].y_g_c)
                         - r(s[0].x_g_n, s[0].y_g_n, s[0].x_g_r, s[0].y_g_r);
            s[0].epsilon_l = s[0].epsilon_nc + MethodData.C * s[0].Delta_r * sign(s[0].dot_epsilon_nc);
            s[0].Theta_r = s[0].epsilon_l + asin(s[0].v_n / s[0].v_r * sin(s[0].Theta_n - s[0].epsilon_l));
            break;

        case GuidanceMethod::Gorka:
            s[0].Theta_r = s[0].epsilon_rc;
            break;

        default:
            break;
        }
    }

    if (isNaN(s[0].Theta_n))
        s[0].Theta_n = s[0].Theta_r;

    // === Начальные перегрузки ===
    s[0].n_xa_c = Linterp(Cdata.n_xa_c_potr, 0.0);
    s[0].n_ya_c = Linterp(Cdata.n_ya_c_potr, 0.0);
    s[0].n_xa_n = Linterp(Ndata.n_xa_n_potr, 0.0);
    s[0].n_ya_n = Linterp(Ndata.n_ya_n_potr, 0.0);
    s[0].n_xa_r = Linterp(Rdata.n_xa_r_potr, 0.0);
    s[0].n_ya_r = 0.0;

    // === Ускорения ===
    double a_xa_c = (s[0].n_xa_c - sin(s[0].Theta_c)) * g;
    double a_ya_c = (s[0].n_ya_c - cos(s[0].Theta_c)) * g;
    double a_xa_n = (s[0].n_xa_n - sin(s[0].Theta_n)) * g;
    double a_ya_n = (s[0].n_ya_n - cos(s[0].Theta_n)) * g;
    double a_xa_r = (s[0].n_xa_r - sin(s[0].Theta_r)) * g;
    double a_ya_r = (s[0].n_ya_r - cos(s[0].Theta_r)) * g;

    // === Цикл интегрирования ===
    unsigned int i = 0;
    std::string returnCode = "t_r_max";

    while (s[i].t < limits.t_r_max && i < N_max - 1)
    {
        if (i >= 1)
        {
            s[i].n_xa_c = Linterp(Cdata.n_xa_c_potr, s[i].t);
            s[i].n_ya_c = Linterp(Cdata.n_ya_c_potr, s[i].t);
            s[i].n_xa_n = Linterp(Ndata.n_xa_n_potr, s[i].t);
            s[i].n_ya_n = Linterp(Ndata.n_ya_n_potr, s[i].t);
            s[i].n_xa_r = Linterp(Rdata.n_xa_r_potr, s[i].t);

            switch (MethodData.Method)
            {
            case GuidanceMethod::Null:
                s[i].n_ya_r = 0.0;
                break;

            case GuidanceMethod::PurePursuit:
            case GuidanceMethod::DeviatedPursuit:
                s[i].n_ya_r = n_ya_potr_PurePursuit(s[i-1].dot_epsilon_rc, s[i-1].v_r, s[i-1].Theta_r);
                break;

            case GuidanceMethod::ProportionalNav_const:
                s[i].n_ya_r = n_ya_potr_Proportional_const(
                    s[i-1].dot_epsilon_rc,
                    s[i-1].v_r,
                    MethodData.k,
                    s[i-1].Theta_r
                );
                break;

            case GuidanceMethod::ProportionalNav_var:
                s[i].n_ya_r = n_ya_potr_Proportional_var(
                    s[i-1].dot_epsilon_rc,
                    s[i-1].dot_r_rc,
                    MethodData.lambda,
                    s[i-1].Theta_r
                );
                break;

            case GuidanceMethod::BeamRiding_LineOfSight:
                s[i-1].epsilon_l = s[i-1].epsilon_nc;
                s[i-1].dot_epsilon_l = s[i-1].dot_epsilon_nc;
                s[i-1].ddot_epsilon_l = s[i-1].ddot_epsilon_nc;
                s[i].n_ya_r = n_ya_potr_BeamRiding(
                    s[i-1].epsilon_l, s[i-1].dot_epsilon_l, s[i-1].ddot_epsilon_l,
                    s[i-1].r_nr, s[i-1].dot_r_nr, s[i-1].Theta_r, s[i-1].Theta_n,
                    a_xa_r, a_xa_n, a_ya_n
                );
                break;

            case GuidanceMethod::BeamRiding_ConstAlignment:
                s[i-1].Delta_r = s[i-1].r_nc - s[i-1].r_nr;
                s[i-1].dot_Delta_r = s[i-1].dot_r_nc - s[i-1].dot_r_nr;
                s[i-1].ddot_Delta_r = s[i-1].ddot_r_nc - s[i-1].ddot_r_nr;
                s[i-1].epsilon_l = s[i-1].epsilon_nc + MethodData.C * sign(s[i-1].dot_epsilon_nc) * s[i-1].Delta_r;
                s[i-1].dot_epsilon_l = s[i-1].dot_epsilon_nc + MethodData.C * sign(s[i-1].dot_epsilon_nc) * s[i-1].dot_Delta_r;
                s[i-1].ddot_epsilon_l = s[i-1].ddot_epsilon_nc + MethodData.C * sign(s[i-1].dot_epsilon_nc) * s[i-1].ddot_Delta_r;
                s[i].n_ya_r = n_ya_potr_BeamRiding(
                    s[i-1].epsilon_l, s[i-1].dot_epsilon_l, s[i-1].ddot_epsilon_l,
                    s[i-1].r_nr, s[i-1].dot_r_nr, s[i-1].Theta_r, s[i-1].Theta_n,
                    a_xa_r, a_xa_n, a_ya_n
                );
                break;

            case GuidanceMethod::EVT:
                s[i].n_ya_r = n_ya_potr_EVT(
                    s[i-1].v_r,
                    s[i-1].dot_epsilon_rc,
                    s[i-1].Theta_r,
                    MethodData.K_g,
                    MethodData.k
                );
                break;

            case GuidanceMethod::Gorka:
            {
                double dTheta_r_prog = 0.0;
                double h_prog_maneuver = 5000.0;
                s[i].n_ya_r = n_ya_potr_Gorka(
                    s[i-1].v_r,
                    s[i-1].Theta_r,
                    s[i-1].y_g_r,
                    s[i-1].r_rc,
                    s[i-1].epsilon_rc,
                    h_prog_maneuver,
                    dTheta_r_prog
                );
                s[i].dTheta_r_prog = dTheta_r_prog;
                break;
            }

            default:
                break;
            }

            a_xa_c = (s[i].n_xa_c - sin(s[i].Theta_c)) * g;
            a_ya_c = (s[i].n_ya_c - cos(s[i].Theta_c)) * g;
            a_xa_n = (s[i].n_xa_n - sin(s[i].Theta_n)) * g;
            a_ya_n = (s[i].n_ya_n - cos(s[i].Theta_n)) * g;
            a_xa_r = (s[i].n_xa_r - sin(s[i].Theta_r)) * g;
            a_ya_r = (s[i].n_ya_r - cos(s[i].Theta_r)) * g;
        }

        // === Кинематика ===
        s[i].r_rc = r(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c);
        s[i].dot_r_rc = dot_r(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c,
                              s[i].v_r, s[i].v_c, s[i].Theta_r, s[i].Theta_c);
        s[i].ddot_r_rc = ddot_r(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c,
                                s[i].v_r, s[i].v_c, s[i].Theta_r, s[i].Theta_c,
                                a_xa_r, a_ya_r, a_xa_c, a_ya_c);
        s[i].epsilon_rc = epsilon(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c);
        s[i].dot_epsilon_rc = dot_epsilon(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c,
                                          s[i].v_r, s[i].v_c, s[i].Theta_r, s[i].Theta_c);
        s[i].ddot_epsilon_rc = ddot_epsilon(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c,
                                            s[i].v_r, s[i].v_c, s[i].Theta_r, s[i].Theta_c,
                                            a_xa_r, a_ya_r, a_xa_c, a_ya_c);

        s[i].r_nr = r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r);
        s[i].dot_r_nr = dot_r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r,
                              s[i].v_n, s[i].v_r, s[i].Theta_n, s[i].Theta_r);
        s[i].ddot_r_nr = ddot_r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r,
                                s[i].v_n, s[i].v_r, s[i].Theta_n, s[i].Theta_r,
                                a_xa_n, a_ya_n, a_xa_r, a_ya_r);
        s[i].epsilon_nr = epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r);
        s[i].dot_epsilon_nr = dot_epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r,
                                          s[i].v_n, s[i].v_r, s[i].Theta_n, s[i].Theta_r);
        s[i].ddot_epsilon_nr = ddot_epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r,
                                            s[i].v_n, s[i].v_r, s[i].Theta_n, s[i].Theta_r,
                                            a_xa_n, a_ya_n, a_xa_r, a_ya_r);

        s[i].r_nc = r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c);
        s[i].dot_r_nc = dot_r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c,
                              s[i].v_n, s[i].v_c, s[i].Theta_n, s[i].Theta_c);
        s[i].ddot_r_nc = ddot_r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c,
                                s[i].v_n, s[i].v_c, s[i].Theta_n, s[i].Theta_c,
                                a_xa_n, a_ya_n, a_xa_c, a_ya_c);
        s[i].epsilon_nc = epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c);
        s[i].dot_epsilon_nc = dot_epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c,
                                          s[i].v_n, s[i].v_c, s[i].Theta_n, s[i].Theta_c);
        s[i].ddot_epsilon_nc = ddot_epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c,
                                            s[i].v_n, s[i].v_c, s[i].Theta_n, s[i].Theta_c,
                                            a_xa_n, a_ya_n, a_xa_c, a_ya_c);

        // === Проверка ограничений ===
        if (fabs(s[i].n_ya_r) > limits.n_ya_r_max)
        {
            returnCode = "n_ya_r_max";
            break;
        }

        if (s[i].t < limits.t_r_min && s[i].r_rc < calcParams.r_por)
        {
            returnCode = "t_r_min";
            break;
        }

        if (s[i].t > limits.t_r_max)
        {
            returnCode = "t_r_max";
            break;
        }

        if (!isNaN(limits.R_r_min) &&
            r(inits.x_g_pusk, inits.y_g_pusk, s[i].x_g_r, s[i].y_g_r) < limits.R_r_min &&
            s[i].r_rc < calcParams.r_por)
        {
            returnCode = "R_r_min";
            break;
        }

        if (!isNaN(limits.R_r_max) &&
            r(inits.x_g_pusk, inits.y_g_pusk, s[i].x_g_r, s[i].y_g_r) > limits.R_r_max)
        {
            returnCode = "R_r_max";
            break;
        }

        if (!isNaN(limits.y_g_r_min) && s[i].y_g_r < limits.y_g_r_min && s[i].t > limits.t_r_min)
        {
            returnCode = "y_g_r_min";
            break;
        }

        if (!isNaN(limits.y_g_r_max) && s[i].y_g_r > limits.y_g_r_max && s[i].t > limits.t_r_min)
        {
            returnCode = "y_g_r_max";
            break;
        }

        if (!isNaN(limits.epsilon_nc_min) &&
            (s[i].epsilon_nc < limits.epsilon_nc_min || s[i].epsilon_nc > M_PI - limits.epsilon_nc_min))
        {
            returnCode = "epsilon_nc_min";
            break;
        }

        if (!isNaN(limits.epsilon_nc_max) &&
            s[i].epsilon_nc > limits.epsilon_nc_max && s[i].epsilon_nc < M_PI - limits.epsilon_nc_max)
        {
            returnCode = "epsilon_nc_max";
            break;
        }

        if (!isNaN(limits.v_r_min) && s[i].v_r < limits.v_r_min)
        {
            returnCode = "v_r_min";
            break;
        }

        if (s[i].r_rc < calcParams.r_por)
        {
            returnCode = "0";
            break;
        }

        if (s[i].y_g_r <= 0.0 && s[i].t > limits.t_r_min)
        {
            returnCode = "ground_hit";
            break;
        }

        // === Шаг Эйлера ===
        TrajectoryParameters next = s[i];
        next.t = s[i].t + dt;

        next.x_g_c += s[i].v_c * cos(s[i].Theta_c) * dt;
        next.y_g_c += s[i].v_c * sin(s[i].Theta_c) * dt;
        next.v_c += a_xa_c * dt;
        if (fabs(s[i].v_c) > 1e-6) next.Theta_c += (a_ya_c / s[i].v_c) * dt;

        next.x_g_n += s[i].v_n * cos(s[i].Theta_n) * dt;
        next.y_g_n += s[i].v_n * sin(s[i].Theta_n) * dt;
        next.v_n += a_xa_n * dt;
        if (fabs(s[i].v_n) > 1e-6) next.Theta_n += (a_ya_n / s[i].v_n) * dt;

        next.x_g_r += s[i].v_r * cos(s[i].Theta_r) * dt;
        next.y_g_r += s[i].v_r * sin(s[i].Theta_r) * dt;
        next.v_r += a_xa_r * dt;
        if (fabs(s[i].v_r) > 1e-6) next.Theta_r += (a_ya_r / s[i].v_r) * dt;

        s.push_back(next);
        ++i;
    }

    return returnCode;
}


// ============================================================================
// 2. ПОЛНАЯ ФИЗИЧЕСКАЯ МОДЕЛЬ (с двигателем, массой, аэродинамикой)
// ============================================================================
std::string CalcTrajectory_DZ(
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
    std::vector<TrajectoryParameters> &s,
    double &m_norm,
    double &switch_time)
{
    double dt = calc.dt;
    unsigned int N_max = static_cast<unsigned int>(ceil(limits.t_r_max / dt));

    // --- Параметры двигателя ---
    const double I_ud = 2450.0;
    double m_pn = rocket.m_pn();
    double diametr = rocket.d;
    double beta = rocket.beta;

    double m_0_start = m_pn / (1.0 - beta * mu_0);
    double m_fuel = mu_0 * m_0_start;

    double m_fuel_raz = m_fuel / (1.0 + K_m);
    double m_fuel_dva = m_fuel_raz * K_m;
    double m_0 = m_0_start;

    double tiaga_raz = eta_0 * m_0 * g;
    double tiaga_dva = tiaga_raz * K_P;

    double G_raz = tiaga_raz / I_ud;
    double G_dva = tiaga_dva / I_ud;

    // --- Инициализация ---
    s.clear();
    s.reserve(N_max);

    TrajectoryParameters p = {};
    p.t = 0.0;
    p.Theta_c = inits.Theta_c0;
    p.v_c = inits.v_c0;
    p.x_g_c = inits.x_g_c0;
    p.y_g_c = inits.y_g_c0;
    p.Theta_n = inits.Theta_n0;
    p.v_n = inits.v_n0;
    p.x_g_n = inits.x_g_pusk;
    p.y_g_n = inits.y_g_pusk;
    p.Theta_r = inits.Theta_r0;
    p.v_r = inits.v_r0;
    p.x_g_r = inits.x_g_pusk;
    p.y_g_r = inits.y_g_pusk;
    p.mass_r = m_0;
    p.P_r = 0.0;
    p.mode1 = true;

    s.push_back(p);

    // === Начальные углы ===
    s[0].epsilon_rc = epsilon(s[0].x_g_r, s[0].y_g_r, s[0].x_g_c, s[0].y_g_c);
    s[0].epsilon_nc = s[0].epsilon_rc;
    s[0].dot_epsilon_rc = dot_epsilon(s[0].x_g_r, s[0].y_g_r, s[0].x_g_c, s[0].y_g_c,
                                      s[0].v_r, s[0].v_c, s[0].epsilon_rc, s[0].Theta_c);
    s[0].dot_epsilon_nc = dot_epsilon(s[0].x_g_n, s[0].y_g_n, s[0].x_g_c, s[0].y_g_c,
                                      s[0].v_n, s[0].v_c, s[0].Theta_n, s[0].Theta_c);

    if (isNaN(s[0].Theta_r))
    {
        switch (MethodData.Method)
        {
        case GuidanceMethod::Null:
        case GuidanceMethod::PurePursuit:
        case GuidanceMethod::ProportionalNav_const:
        case GuidanceMethod::ProportionalNav_var:
        case GuidanceMethod::EVT:
            s[0].Theta_r = s[0].epsilon_rc;
            break;

        case GuidanceMethod::DeviatedPursuit:
            s[0].Theta_r = s[0].epsilon_rc + MethodData.phi_upr * sign(s[0].dot_epsilon_rc);
            break;

        case GuidanceMethod::ParallelNav:
            s[0].Theta_r = s[0].epsilon_rc - asin(s[0].v_c / s[0].v_r * sin(s[0].epsilon_rc - s[0].Theta_c));
            break;

        case GuidanceMethod::BeamRiding_LineOfSight:
            s[0].epsilon_l = s[0].epsilon_nc;
            s[0].Theta_r = s[0].epsilon_l + asin(s[0].v_n / s[0].v_r * sin(s[0].Theta_n - s[0].epsilon_l));
            break;

        case GuidanceMethod::BeamRiding_ConstAlignment:
            s[0].Delta_r = r(s[0].x_g_n, s[0].y_g_n, s[0].x_g_c, s[0].y_g_c)
                         - r(s[0].x_g_n, s[0].y_g_n, s[0].x_g_r, s[0].y_g_r);
            s[0].epsilon_l = s[0].epsilon_nc + MethodData.C * s[0].Delta_r * sign(s[0].dot_epsilon_nc);
            s[0].Theta_r = s[0].epsilon_l + asin(s[0].v_n / s[0].v_r * sin(s[0].Theta_n - s[0].epsilon_l));
            break;

        case GuidanceMethod::Gorka:
            s[0].Theta_r = s[0].epsilon_rc;
            break;

        default:
            break;
        }
    }

    if (isNaN(s[0].Theta_n))
        s[0].Theta_n = s[0].Theta_r;

    // === Начальные перегрузки ===
    s[0].n_xa_c = Linterp(Cdata.n_xa_c_potr, 0.0);
    s[0].n_ya_c = Linterp(Cdata.n_ya_c_potr, 0.0);
    s[0].n_xa_n = Linterp(Ndata.n_xa_n_potr, 0.0);
    s[0].n_ya_n = Linterp(Ndata.n_ya_n_potr, 0.0);
    s[0].n_xa_r = Linterp(Rdata.n_xa_r_potr, 0.0);
    s[0].n_ya_r = 0.0;

    // === Ускорения ===
    double a_xa_c = (s[0].n_xa_c - sin(s[0].Theta_c)) * g;
    double a_ya_c = (s[0].n_ya_c - cos(s[0].Theta_c)) * g;
    double a_xa_n = (s[0].n_xa_n - sin(s[0].Theta_n)) * g;
    double a_ya_n = (s[0].n_ya_n - cos(s[0].Theta_n)) * g;
    double a_xa_r = (s[0].n_xa_r - sin(s[0].Theta_r)) * g;
    double a_ya_r = (s[0].n_ya_r - cos(s[0].Theta_r)) * g;

    // === Цикл интегрирования ===
    unsigned int i = 0;
    std::string returnCode = "t_r_max";

    while (s[i].t < limits.t_r_max && i < N_max - 1)
    {
        // --- 1. Атмосфера ---
        double rho = Atm_GOST4401::rho_atm(s[i].y_g_r);
        double a_snd = Atm_GOST4401::a_atm(s[i].y_g_r);
        double M = (a_snd > 1e-6) ? (s[i].v_r / a_snd) : 0.0;

        // --- 2. Двигатель ---
        double P_val = 0.0;
        bool currentMode1 = false;

        if (m_fuel >= m_fuel_dva && m_fuel_raz != 0.0)
        {
            P_val = tiaga_raz + (Atm_GOST4401::p_atm(0.0) - Atm_GOST4401::p_atm(s[i].y_g_r)) * S_a(diametr);
            currentMode1 = true;
        }
        else if (m_fuel > 0.0)
        {
            P_val = tiaga_dva + (Atm_GOST4401::p_atm(0.0) - Atm_GOST4401::p_atm(s[i].y_g_r)) * S_a(diametr);
            currentMode1 = false;
        }

        // --- 3. Расход топлива ---
        if (m_fuel >= m_fuel_dva)
            m_fuel -= G_raz * dt;
        else if (m_fuel > 0.0)
            m_fuel -= G_dva * dt;

        if (m_fuel <= 0.0)
        {
            switch_time = s[i].t;
            m_fuel = -1.0;
        }

        s[i].mass_r = m_0;
        s[i].P_r = P_val;
        s[i].mode1 = currentMode1;

        // --- 4. Кинематика ---
        s[i].r_rc = r(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c);
        s[i].epsilon_rc = epsilon(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c);

        if (s[i].r_rc > 1e-6)
            s[i].dot_epsilon_rc = dot_epsilon(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c,
                                              s[i].v_r, s[i].v_c, s[i].Theta_r, s[i].Theta_c);
        else
            s[i].dot_epsilon_rc = 0.0;

        s[i].r_nr = r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r);
        s[i].dot_r_nr = dot_r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r,
                              s[i].v_n, s[i].v_r, s[i].Theta_n, s[i].Theta_r);
        s[i].epsilon_nr = epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r);
        s[i].dot_epsilon_nr = dot_epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r,
                                          s[i].v_n, s[i].v_r, s[i].Theta_n, s[i].Theta_r);

        s[i].r_nc = r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c);
        s[i].dot_r_nc = dot_r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c,
                              s[i].v_n, s[i].v_c, s[i].Theta_n, s[i].Theta_c);
        s[i].epsilon_nc = epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c);
        s[i].dot_epsilon_nc = dot_epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c,
                                          s[i].v_n, s[i].v_c, s[i].Theta_n, s[i].Theta_c);

        // --- 5. Наведение ---
        if (i >= 1)
        {
            switch (MethodData.Method)
            {
            case GuidanceMethod::Null:
                s[i].n_ya_r = 0.0;
                break;

            case GuidanceMethod::PurePursuit:
            case GuidanceMethod::DeviatedPursuit:
                s[i].n_ya_r = n_ya_potr_PurePursuit(s[i-1].dot_epsilon_rc, s[i-1].v_r, s[i-1].Theta_r);
                break;

            case GuidanceMethod::ProportionalNav_const:
                s[i].n_ya_r = n_ya_potr_Proportional_const(
                    s[i-1].dot_epsilon_rc,
                    s[i-1].v_r,
                    MethodData.k,
                    s[i-1].Theta_r
                );
                break;

            case GuidanceMethod::ProportionalNav_var:
                s[i].n_ya_r = n_ya_potr_Proportional_var(
                    s[i-1].dot_epsilon_rc,
                    s[i-1].dot_r_rc,
                    MethodData.lambda,
                    s[i-1].Theta_r
                );
                break;

            case GuidanceMethod::BeamRiding_LineOfSight:
                s[i-1].epsilon_l = s[i-1].epsilon_nc;
                s[i-1].dot_epsilon_l = s[i-1].dot_epsilon_nc;
                s[i-1].ddot_epsilon_l = s[i-1].ddot_epsilon_nc;
                s[i].n_ya_r = n_ya_potr_BeamRiding(
                    s[i-1].epsilon_l, s[i-1].dot_epsilon_l, s[i-1].ddot_epsilon_l,
                    s[i-1].r_nr, s[i-1].dot_r_nr, s[i-1].Theta_r, s[i-1].Theta_n,
                    a_xa_r, a_xa_n, a_ya_n
                );
                break;

            case GuidanceMethod::BeamRiding_ConstAlignment:
                s[i-1].Delta_r = s[i-1].r_nc - s[i-1].r_nr;
                s[i-1].dot_Delta_r = s[i-1].dot_r_nc - s[i-1].dot_r_nr;
                s[i-1].ddot_Delta_r = s[i-1].ddot_r_nc - s[i-1].ddot_r_nr;
                s[i-1].epsilon_l = s[i-1].epsilon_nc + MethodData.C * sign(s[i-1].dot_epsilon_nc) * s[i-1].Delta_r;
                s[i-1].dot_epsilon_l = s[i-1].dot_epsilon_nc + MethodData.C * sign(s[i-1].dot_epsilon_nc) * s[i-1].dot_Delta_r;
                s[i-1].ddot_epsilon_l = s[i-1].ddot_epsilon_nc + MethodData.C * sign(s[i-1].dot_epsilon_nc) * s[i-1].ddot_Delta_r;
                s[i].n_ya_r = n_ya_potr_BeamRiding(
                    s[i-1].epsilon_l, s[i-1].dot_epsilon_l, s[i-1].ddot_epsilon_l,
                    s[i-1].r_nr, s[i-1].dot_r_nr, s[i-1].Theta_r, s[i-1].Theta_n,
                    a_xa_r, a_xa_n, a_ya_n
                );
                break;

            case GuidanceMethod::EVT:
                s[i].n_ya_r = n_ya_potr_EVT(
                    s[i-1].v_r,
                    s[i-1].dot_epsilon_rc,
                    s[i-1].Theta_r,
                    MethodData.K_g,
                    MethodData.k
                );
                break;

            case GuidanceMethod::Gorka:
            {
                double dTheta_r_prog = 0.0;
                double h_prog_maneuver = 5000.0;
                s[i].n_ya_r = n_ya_potr_Gorka(
                    s[i-1].v_r,
                    s[i-1].Theta_r,
                    s[i-1].y_g_r,
                    s[i-1].r_rc,
                    s[i-1].epsilon_rc,
                    h_prog_maneuver,
                    dTheta_r_prog
                );
                s[i].dTheta_r_prog = dTheta_r_prog;
                break;
            }

            default:
                break;
            }
        }

        // --- 6. Аэродинамика ---
        double alpha_r = 0.0;
        double delta_r = 0.0;

        if (s[i].v_r > 1.0)
        {
            alpha_r = find_alpha(M, rho, s[i].v_r, diametr, m_0, P_val, s[i].n_ya_r);
            delta_r = delta_bal_over_alpha(M) * alpha_r;
        }

        s[i].alpha_r = alpha_r;
        s[i].delta_r = delta_r;

        double Xa_val = X_a(M, alpha_r, delta_r, rho, s[i].v_r, diametr, P_val);
        s[i].n_xa_r = (P_val * std::cos(alpha_r) - Xa_val) / (m_0 * g);

        // --- 7. Проверки ---
        if (std::fabs(alpha_r) >= 20.0 * M_PI / 180.0 - 1e-6)
        {
            returnCode = "alpha_max";
            break;
        }

        if (fabs(s[i].n_ya_r) > limits.n_ya_r_max)
        {
            returnCode = "n_ya_r_max";
            break;
        }

        if (fabs(s[i].n_xa_r) > limits.n_ya_r_max)
        {
            returnCode = "n_xa_r_max";
            break;
        }

        if (s[i].t < limits.t_r_min && s[i].r_rc < calc.r_por)
        {
            returnCode = "t_r_min";
            break;
        }

        if (s[i].t > limits.t_r_max)
        {
            returnCode = "t_r_max";
            break;
        }

        if (!isNaN(limits.R_r_min) &&
            r(inits.x_g_pusk, inits.y_g_pusk, s[i].x_g_r, s[i].y_g_r) < limits.R_r_min &&
            s[i].r_rc < calc.r_por)
        {
            returnCode = "R_r_min";
            break;
        }

        if (!isNaN(limits.R_r_max) &&
            r(inits.x_g_pusk, inits.y_g_pusk, s[i].x_g_r, s[i].y_g_r) > limits.R_r_max)
        {
            returnCode = "R_r_max";
            break;
        }

        if (!isNaN(limits.y_g_r_min) && s[i].y_g_r < limits.y_g_r_min && s[i].t > limits.t_r_min)
        {
            returnCode = "y_g_r_min";
            break;
        }

        if (!isNaN(limits.y_g_r_max) && s[i].y_g_r > limits.y_g_r_max && s[i].t > limits.t_r_min)
        {
            returnCode = "y_g_r_max";
            break;
        }

        if (!isNaN(limits.epsilon_nc_min) &&
            (s[i].epsilon_nc < limits.epsilon_nc_min || s[i].epsilon_nc > M_PI - limits.epsilon_nc_min))
        {
            returnCode = "epsilon_nc_min";
            break;
        }

        if (!isNaN(limits.epsilon_nc_max) &&
            s[i].epsilon_nc > limits.epsilon_nc_max && s[i].epsilon_nc < M_PI - limits.epsilon_nc_max)
        {
            returnCode = "epsilon_nc_max";
            break;
        }

        if (!isNaN(limits.v_r_min) && s[i].v_r < limits.v_r_min)
        {
            returnCode = "v_r_min";
            break;
        }

        if (s[i].r_rc < calc.r_por)
        {
            returnCode = "0";
            break;
        }

        if (s[i].y_g_r <= 0.0 && s[i].t > limits.t_r_min)
        {
            returnCode = "ground_hit";
            break;
        }

        // --- 8. Ускорения ---
        a_xa_r = (s[i].n_xa_r - sin(s[i].Theta_r)) * g;
        a_ya_r = (s[i].n_ya_r - cos(s[i].Theta_r)) * g;

        // --- 9. Шаг Эйлера ---
        TrajectoryParameters next = s[i];
        next.t = s[i].t + dt;
        next.x_g_r += s[i].v_r * cos(s[i].Theta_r) * dt;
        next.y_g_r += s[i].v_r * sin(s[i].Theta_r) * dt;
        next.v_r += a_xa_r * dt;
        if (fabs(s[i].v_r) > 1e-6) next.Theta_r += (a_ya_r / s[i].v_r) * dt;

        // Цель
        next.x_g_c += s[i].v_c * cos(s[i].Theta_c) * dt;
        next.y_g_c += s[i].v_c * sin(s[i].Theta_c) * dt;
        next.v_c += (s[i].n_xa_c - sin(s[i].Theta_c)) * g * dt;
        if (fabs(s[i].v_c) > 1e-6) next.Theta_c += ((s[i].n_ya_c - cos(s[i].Theta_c)) * g / s[i].v_c) * dt;

        // Носитель
        next.x_g_n += s[i].v_n * cos(s[i].Theta_n) * dt;
        next.y_g_n += s[i].v_n * sin(s[i].Theta_n) * dt;
        next.v_n += (s[i].n_xa_n - sin(s[i].Theta_n)) * g * dt;
        if (fabs(s[i].v_n) > 1e-6) next.Theta_n += ((s[i].n_ya_n - cos(s[i].Theta_n)) * g / s[i].v_n) * dt;

        s.push_back(next);
        ++i;
    }

    // --- Пост-обработка ---
    if (s.size() >= 2)
    {
        s[0].mass_r = s[1].mass_r;
        s[0].P_r = s[1].P_r;
    }

    m_norm = m_fuel / m_0_start;

    return returnCode;
}