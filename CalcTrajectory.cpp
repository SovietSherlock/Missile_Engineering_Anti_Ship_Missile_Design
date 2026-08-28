#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <limits>
#include "CalcTrajectory.h"
#include "DataStructs.h"
#include "LibConstFunc.h"
#include "KinematicFunc2D.h"
#include "Atm_GOST4401.h"
#include "GuidanceMethods.h"

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
    std::vector<TrajectoryParameters> &s)
{
    double dt = calc.dt;
    unsigned int N_max = static_cast<unsigned int>(std::ceil(limits.t_r_max / dt));

    // --- Параметры двигателя ---
    const double I_ud = 2450.0;
    double m_pn    = rocket.m_pn();     // масса полезной нагрузки (БЧ + ОУН + НП)
    double diametr = rocket.d;
    double beta    = rocket.beta;

    // m_0 = m_пн / (1 - beta·mu_0),  m_t = mu_0·m_0
    double m_0_start = m_pn / (1.0 - beta * mu_0);
    double m_fuel    = mu_0 * m_0_start;        // стартовая масса топлива
    double m_sukh    = m_0_start - m_fuel;      // сухая масса (конструкция + БЧ + ОУН + НП)

    double m_fuel_raz = m_fuel / (1.0 + K_m);
    double m_fuel_dva = m_fuel_raz * K_m;
    double m_0 = m_0_start;

    double tiaga_raz = eta_0 * m_0 * g;
    double tiaga_dva = tiaga_raz * K_P;

    double G_raz = tiaga_raz / I_ud;
    double G_dva = tiaga_dva / I_ud;

    // --- Инициализация траектории ---
    s.clear();
    s.reserve(static_cast<int>(N_max));

    TrajectoryParameters p = {};

    p.t = 0.0;
    p.Theta_c = inits.Theta_c0;  p.v_c = inits.v_c0;
    p.x_g_c = inits.x_g_c0;      p.y_g_c = inits.y_g_c0;
    p.Theta_n = inits.Theta_n0;  p.v_n = inits.v_n0;
    p.x_g_n = inits.x_g_pusk;    p.y_g_n = inits.y_g_pusk;
    p.Theta_r = inits.Theta_r0;  p.v_r = inits.v_r0;
    p.x_g_r = inits.x_g_pusk;    p.y_g_r = inits.y_g_pusk;

    p.mass_r = m_0;
    p.P_r = 0.0;
    p.mode1 = true;

    s.push_back(p);

    // === Начальные углы ===

    s[0].epsilon_rc = epsilon(s[0].x_g_r, s[0].y_g_r, s[0].x_g_c, s[0].y_g_c);
    s[0].epsilon_nc = s[0].epsilon_rc;

    if (std::isnan(s[0].Theta_r))
    {
        switch (MethodData.Method)
        {
        case GuidanceMethod::Null:
        case GuidanceMethod::PurePursuit:
        case GuidanceMethod::ProportionalNav_const:
        case GuidanceMethod::ProportionalNav_var:
        case GuidanceMethod::EVT:
        case GuidanceMethod::Gorka:
            s[0].Theta_r = s[0].epsilon_rc;
            break;

        case GuidanceMethod::DeviatedPursuit:
            s[0].dot_epsilon_rc = dot_epsilon(s[0].x_g_r, s[0].y_g_r, s[0].x_g_c, s[0].y_g_c,
                                              s[0].v_r, s[0].v_c, s[0].epsilon_rc, s[0].Theta_c);
            s[0].Theta_r = s[0].epsilon_rc + MethodData.phi_upr * sign(s[0].dot_epsilon_rc);
            break;

        case GuidanceMethod::ParallelNav:
            s[0].Theta_r = s[0].epsilon_rc - std::asin(s[0].v_c / s[0].v_r * std::sin(s[0].epsilon_rc - s[0].Theta_c));
            break;

        case GuidanceMethod::BeamRiding_LineOfSight:
            s[0].epsilon_l = s[0].epsilon_nc;
            s[0].Theta_r = s[0].epsilon_l + std::asin(s[0].v_n / s[0].v_r * std::sin(s[0].Theta_n - s[0].epsilon_l));
            break;

        case GuidanceMethod::BeamRiding_ConstAlignment:
            s[0].Delta_r = r(s[0].x_g_n, s[0].y_g_n, s[0].x_g_c, s[0].y_g_c)
                           - r(s[0].x_g_n, s[0].y_g_n, s[0].x_g_r, s[0].y_g_r);
            s[0].dot_epsilon_nc = dot_epsilon(s[0].x_g_n, s[0].y_g_n, s[0].x_g_c, s[0].y_g_c,
                                              s[0].v_n, s[0].v_c, s[0].Theta_n, s[0].Theta_c);
            s[0].epsilon_l = s[0].epsilon_nc + MethodData.C * s[0].Delta_r * sign(s[0].dot_epsilon_nc);
            s[0].Theta_r = s[0].epsilon_l + std::asin(s[0].v_n / s[0].v_r * std::sin(s[0].Theta_n - s[0].epsilon_l));
            break;

        default:
            break;
        }
    }

    if (std::isnan(s[0].Theta_n))
        s[0].Theta_n = s[0].Theta_r;

    // === Начальные перегрузки ===
    s[0].n_xa_c = Linterp(Cdata.n_xa_c_potr, 0.0);
    s[0].n_ya_c = Linterp(Cdata.n_ya_c_potr, 0.0);
    s[0].n_xa_n = Linterp(Ndata.n_xa_n_potr, 0.0);
    s[0].n_ya_n = Linterp(Ndata.n_ya_n_potr, 0.0);
    s[0].n_xa_r = Linterp(Rdata.n_xa_r_potr, 0.0);
    s[0].n_ya_r = 0.0;

    // === Ускорения ===
    double a_xa_c = (s[0].n_xa_c - std::sin(s[0].Theta_c)) * g;
    double a_ya_c = (s[0].n_ya_c - std::cos(s[0].Theta_c)) * g;
    double a_xa_n = (s[0].n_xa_n - std::sin(s[0].Theta_n)) * g;
    double a_ya_n = (s[0].n_ya_n - std::cos(s[0].Theta_n)) * g;
    double a_xa_r = (s[0].n_xa_r - std::sin(s[0].Theta_r)) * g;
    double a_ya_r = (s[0].n_ya_r - std::cos(s[0].Theta_r)) * g;

    // === Цикл интегрирования ===
    unsigned int i = 0;
    std::string returnCode = "t_r_max";

    // === Номер фазы ===
    int phase = 0; // Начальная фаза для многофазного метода (0 = ЭВТ на виртуальную цель, 1 = Маршевый полет, 2 = Полет по "горке", 3 = Пикирование)

    // === Критическое расстояние для пикирования ===
    double r_crit_calc = calc_r_critical(s[i].v_r, s[i].v_c, MethodData.H_gorka, MethodData.k_dive, limits.n_ya_r_max);

    // === Флаг для отсечения выныривания радиогоризонта ===
    bool radar_horizon_ended = false;

    // === Флаг отсечения учета колебаний траектории на маршевой высоте ==
    bool reached_march = false;


    while (s[i].t < limits.t_r_max && i < N_max - 1)
    {

        // --- 1. Атмосфера ---
        double rho   = Atm_GOST4401::rho_atm(s[i].y_g_r);
        double a_snd = Atm_GOST4401::a_atm(s[i].y_g_r);
        double M     = (a_snd > 1e-6) ? (s[i].v_r / a_snd) : 0.0;

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
            m_fuel = -1.0;

        }

        s[i].mass_r = m_0;
        s[i].P_r    = P_val;
        s[i].mode1  = currentMode1;

        // --- 4. Кинематика (для наведения) ---
        //Расстояние между ракетой и целью, м
        s[i].r_rc = r(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c);
        //критическое расстояние между ракетой и целью для начала пикирования
        double r_crit_calc = calc_r_critical(s[i].v_r, s[i].v_c, MethodData.H_gorka, MethodData.k_dive, limits.n_ya_r_max);
        //Скорость изменения расстояния между ракетой и целью, м/с
        s[i].dot_r_rc = dot_r(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c,
                              s[i].v_r, s[i].v_c, s[i].Theta_r, s[i].Theta_c);
        //Ускорение изменения расстояния между ракетой и целью, м/с^2
        s[i].ddot_r_rc = ddot_r(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c,
                                s[i].v_r, s[i].v_c, s[i].Theta_r, s[i].Theta_c,
                                a_xa_r, a_ya_r, a_xa_c, a_ya_c);
        //Угол наклона линии ракета-цель, рад
        s[i].epsilon_rc = epsilon(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c, s[i].y_g_c);
        //Угловая скорость линии ракета-цель, рад/с
        s[i].dot_epsilon_rc = dot_epsilon(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c,
                                          s[i].y_g_c, s[i].v_r, s[i].v_c, s[i].Theta_r, s[i].Theta_c);
        //угловая скорость линии визирования ракета-виртуальная цель, рад/с
        s[i].dot_epsilon_virt = dot_epsilon_virtual(
            s[i].x_g_r, s[i].y_g_r,
            s[i].x_g_c, s[i].y_g_c,
            MethodData.delta_x_virt, 
            s[i].v_r, s[i].v_c,
            s[i].Theta_r, s[i].Theta_c
        );
        //Угловое ускорение линии ракета-цель, рад/с^2
        s[i].ddot_epsilon_rc = ddot_epsilon(s[i].x_g_r, s[i].y_g_r, s[i].x_g_c,
                                            s[i].y_g_c, s[i].v_r, s[i].v_c, s[i].Theta_r, s[i].Theta_c,
                                            a_xa_r, a_ya_r, a_xa_c, a_ya_c);

        //Расстояние между носителем и ракетой, м
        s[i].r_nr = r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r);
        //Скорость изменения расстояния между носителем и ракетой, м/с
        s[i].dot_r_nr = dot_r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r,
                              s[i].v_n, s[i].v_r, s[i].Theta_n, s[i].Theta_r);
        //Ускорение изменения расстояния между носителем и ракетой, м/с^2
        s[i].ddot_r_nr = ddot_r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r,
                                s[i].v_n, s[i].v_r, s[i].Theta_n, s[i].Theta_r,
                                a_xa_n, a_ya_n, a_xa_r, a_ya_r);
        //Угол наклона линии носитель-ракета, рад
        s[i].epsilon_nr = epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r, s[i].y_g_r);
        //Угловая скорость линии носитель-ракета, рад/с
        s[i].dot_epsilon_nr = dot_epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r,
                                          s[i].y_g_r, s[i].v_n, s[i].v_r, s[i].Theta_n, s[i].Theta_r);
        //Угловое ускорение линии носитель-ракета, рад/с^2
        s[i].ddot_epsilon_nr = ddot_epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_r,
                                            s[i].y_g_r, s[i].v_n, s[i].v_r, s[i].Theta_n, s[i].Theta_r,
                                            a_xa_n, a_ya_n, a_xa_r, a_ya_r);

        //Расстояние между носителем и целью, м
        s[i].r_nc = r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c);
        //Скорость изменения расстояния между носителем и целью, м/с
        s[i].dot_r_nc = dot_r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c,
                              s[i].v_n, s[i].v_c, s[i].Theta_n, s[i].Theta_c);
        //Ускорение изменения расстояния между носителем и целью, м/с^2
        s[i].ddot_r_nc = ddot_r(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c,
                                s[i].v_n, s[i].v_c, s[i].Theta_n, s[i].Theta_c,
                                a_xa_n, a_ya_n, a_xa_c, a_ya_c);
        //Угол наклона линии носитель-цель, рад
        s[i].epsilon_nc = epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c, s[i].y_g_c);
        //Угловая скорость линии носитель-цель, рад/с
        s[i].dot_epsilon_nc = dot_epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c,
                                          s[i].y_g_c, s[i].v_n, s[i].v_c, s[i].Theta_n, s[i].Theta_c);
        //Угловое ускорение линии носитель-цель, рад/с^2
        s[i].ddot_epsilon_nc = ddot_epsilon(s[i].x_g_n, s[i].y_g_n, s[i].x_g_c,
                                            s[i].y_g_c, s[i].v_n, s[i].v_c, s[i].Theta_n, s[i].Theta_c,
                                            a_xa_n, a_ya_n, a_xa_c, a_ya_c);
        //Критическое расстояние пикирования для многофазной траетории, м
        r_crit_calc = calc_r_critical(
                    s[i].v_r, s[i].v_c,
                    MethodData.H_gorka,
                    MethodData.k_dive,
                    limits.n_ya_r_max);

        double y_bez_temp = calc_radar_horizon_boundary(
            s[i].r_rc, 
            MethodData.H_ant,    
            MethodData.H_save,     
            MethodData.H_march
        );

        // Если флаг уже true, значит мы уже под землей. Оставляем nullopt.
        if (radar_horizon_ended)
        {
            s[i].y_radar_horizon = std::nullopt;
        }
        else
        {
            if (y_bez_temp <= 0.0)
            {
                radar_horizon_ended = true; // Запоминаем момент ухода под землю
                s[i].y_radar_horizon = std::nullopt;
            }
            else
            {
                s[i].y_radar_horizon = y_bez_temp;
            }
        }


        // --- 5. Наведение ---

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
                s[i].n_ya_r = n_ya_potr_PurePursuit(s[i-1].dot_epsilon_rc, s[i-1].v_r, s[i-1].Theta_r);
                break;

            case GuidanceMethod::DeviatedPursuit:
                s[i].n_ya_r = n_ya_potr_PurePursuit(s[i-1].dot_epsilon_rc, s[i-1].v_r, s[i-1].Theta_r);
                break;

            case GuidanceMethod::BeamRiding_LineOfSight:
                s[i-1].epsilon_l = s[i-1].epsilon_nc;
                s[i-1].dot_epsilon_l = s[i-1].dot_epsilon_nc;
                s[i-1].ddot_epsilon_l = s[i-1].ddot_epsilon_nc;
                s[i].n_ya_r = n_ya_potr_BeamRiding(s[i-1].epsilon_l, s[i-1].dot_epsilon_l, s[i-1].ddot_epsilon_l,
                                                   s[i-1].r_nr, s[i-1].dot_r_nr, s[i-1].Theta_r, s[i-1].Theta_n,
                                                   a_xa_r, a_xa_n, a_ya_n);
                break;

            case GuidanceMethod::BeamRiding_ConstAlignment:
                s[i-1].Delta_r = s[i-1].r_nc - s[i-1].r_nr;
                s[i-1].dot_Delta_r = s[i-1].dot_r_nc - s[i-1].dot_r_nr;
                s[i-1].ddot_Delta_r = s[i-1].ddot_r_nc - s[i-1].ddot_r_nr;
                s[i-1].epsilon_l = s[i-1].epsilon_nc + MethodData.C * sign(s[i-1].dot_epsilon_nc) * s[i-1].Delta_r;
                s[i-1].dot_epsilon_l = s[i-1].dot_epsilon_nc + MethodData.C * sign(s[i-1].dot_epsilon_nc) * s[i-1].dot_Delta_r;
                s[i-1].ddot_epsilon_l = s[i-1].ddot_epsilon_nc + MethodData.C * sign(s[i-1].dot_epsilon_nc) * s[i-1].ddot_Delta_r;
                s[i].n_ya_r = n_ya_potr_BeamRiding(s[i-1].epsilon_l, s[i-1].dot_epsilon_l, s[i-1].ddot_epsilon_l,
                                                   s[i-1].r_nr, s[i-1].dot_r_nr, s[i-1].Theta_r, s[i-1].Theta_n,
                                                   a_xa_r, a_xa_n, a_ya_n);
                break;

            case GuidanceMethod::ProportionalNav_var:
                s[i].n_ya_r = n_ya_potr_Proportional(MethodData.k, s[i-1].dot_epsilon_rc,
                                                     s[i-1].v_r,
                                                     s[i-1].Theta_r);
                break;

            case GuidanceMethod::ProportionalNav_const:
                s[i].n_ya_r = n_ya_potr_Proportional(MethodData.k, s[i-1].dot_epsilon_rc,
                                                     s[i-1].v_r,
                                                     s[i-1].Theta_r);
                break;

            case GuidanceMethod::EVT:
                s[i].n_ya_r = n_ya_potr_EVT(s[i].v_r, s[i].dot_epsilon_rc, s[i].Theta_r, K_g, MethodData.k);
                break;

            case GuidanceMethod::Gorka:
            {
                // ВЫЗОВ МНОГОФАЗНОЙ ФУНКЦИИ НАВЕДЕНИЯ
                s[i].n_ya_r = n_ya_potr_TwoPhase(
                    s[i].v_r, s[i].Theta_r,
                    s[i].y_g_r,
                    s[i].r_rc,
                    s[i].epsilon_rc,
                    s[i].v_c, s[i].Theta_c,
                    s[i].dot_epsilon_virt, // Угловая скорость на виртуальную цель (посчитана выше)
                    s[i].dot_epsilon_rc,      // Угловая скорость на реальную цель
                    MethodData,
                    phase,                    // Переменная фазы (передается по ссылке, будет меняться внутри)              // Критическое расстояние для пикирования (посчитано выше)
                    limits.n_ya_r_max         // Максимальная допустимая перегрузка
                );

                // std::cout << ", r_crit_calc=" << r_crit_calc;
                // std::cout << ", r_rc=" << s[i].r_rc;
                break;
            }


            default:
                break;
            }


            a_xa_c = (s[i].n_xa_c - std::sin(s[i].Theta_c)) * g;
            a_ya_c = (s[i].n_ya_c - std::cos(s[i].Theta_c)) * g;
            a_xa_n = (s[i].n_xa_n - std::sin(s[i].Theta_n)) * g;
            a_ya_n = (s[i].n_ya_n - std::cos(s[i].Theta_n)) * g;
            a_xa_r = (s[i].n_xa_r - std::sin(s[i].Theta_r)) * g;
            a_ya_r = (s[i].n_ya_r - std::cos(s[i].Theta_r)) * g;

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
        if (std::abs(s[i].n_ya_r) > limits.n_ya_r_max)
        {
            returnCode = "n_ya_r_max";
            break;
        }

        if (s[i].t < limits.t_r_min && s[i].r_rc < calc.r_por)
        {
            returnCode = "t_r_min";
            break;
        }

        if (!std::isnan(limits.t_r_max) && s[i].t > limits.t_r_max)
        {
            returnCode = "t_r_max";
            break;
        }

        if (!std::isnan(limits.R_r_min) &&
            r(inits.x_g_pusk, inits.y_g_pusk, s[i].x_g_r, s[i].y_g_r) < limits.R_r_min &&
            s[i].r_rc < calc.r_por)
        {
            returnCode = "R_r_min";
            break;
        }

        if (!std::isnan(limits.R_r_max) &&
            r(inits.x_g_pusk, inits.y_g_pusk, s[i].x_g_r, s[i].y_g_r) > limits.R_r_max)
        {
            returnCode = "R_r_max";
            break;
        }

        if (!std::isnan(limits.y_g_r_min) && s[i].y_g_r < limits.y_g_r_min && s[i].t > limits.t_r_min)
        {
            returnCode = "y_g_r_min";
            break;
        }

        if (!std::isnan(limits.y_g_r_max) && s[i].y_g_r > limits.y_g_r_max && s[i].t > limits.t_r_min)
        {
            returnCode = "y_g_r_max";
            break;
        }

        if (!std::isnan(limits.epsilon_nc_min) &&
            (s[i].epsilon_nc < limits.epsilon_nc_min ||
             s[i].epsilon_nc > M_PI - limits.epsilon_nc_min))
        {
            returnCode = "epsilon_nc_min";
            break;
        }

        if (!std::isnan(limits.epsilon_nc_max) &&
            s[i].epsilon_nc > limits.epsilon_nc_max &&
            s[i].epsilon_nc < M_PI - limits.epsilon_nc_max)
        {
            returnCode = "epsilon_nc_max";
            break;
        }

        if (s[i].v_r < 119.0)
        {
            returnCode = "v_r_min";
            break;
        }

        // // Проверка радиотени (только до первого касания маршевой высоты)
        // if (MethodData.Method == GuidanceMethod::Gorka && !reached_march && 
        //     s[i].y_radar_horizon != std::nullopt && 
        //     s[i].y_g_r > s[i].y_radar_horizon)
        // {
        //     returnCode = "radar_visibility";
        //     break;
        // }


        if (s[i].r_rc < calc.r_por)
        {
            if (s[i].v_r < 280.0)
            {
                returnCode = "v_r_destruct";
                break;
            }
            else
            {
                returnCode = "0";
                break;
            }
        }

        // --- 2.4. Вычисление производных (правых частей системы диф.ур.) ---
        double dot_Theta_c = 0.0;
        if (s[i].v_c != 0.0)
            dot_Theta_c = (s[i].n_ya_c - std::cos(s[i].Theta_c)) * g / s[i].v_c;
        double dot_v_c = (s[i].n_xa_c - std::sin(s[i].Theta_c)) * g;
        double dot_x_g_c = s[i].v_c * std::cos(s[i].Theta_c);
        double dot_y_g_c = s[i].v_c * std::sin(s[i].Theta_c);

        double dot_Theta_n = 0.0;
        if (s[i].v_n != 0.0)
            dot_Theta_n = (s[i].n_ya_n - std::cos(s[i].Theta_n)) * g / s[i].v_n;
        double dot_v_n = (s[i].n_xa_n - std::sin(s[i].Theta_n)) * g;
        double dot_x_g_n = s[i].v_n * std::cos(s[i].Theta_n);
        double dot_y_g_n = s[i].v_n * std::sin(s[i].Theta_n);

        double dot_Theta_r = 0.0;
        if (s[i].v_r != 0.0)
            dot_Theta_r = (s[i].n_ya_r - std::cos(s[i].Theta_r)) * g / s[i].v_r;
        double dot_v_r = (s[i].n_xa_r - std::sin(s[i].Theta_r)) * g;
        double dot_x_g_r = s[i].v_r * std::cos(s[i].Theta_r);
        double dot_y_g_r = s[i].v_r * std::sin(s[i].Theta_r);

        // --- 2.5. Вычисление значений переменных системы диф.ур. на следующем шаге ---
        TrajectoryParameters next = s[i];

        next.Theta_c = s[i].Theta_c + dot_Theta_c * dt;
        next.v_c     = s[i].v_c     + dot_v_c     * dt;
        next.x_g_c   = s[i].x_g_c   + dot_x_g_c   * dt;
        next.y_g_c   = s[i].y_g_c   + dot_y_g_c   * dt;

        next.Theta_n = s[i].Theta_n + dot_Theta_n * dt;
        next.v_n     = s[i].v_n     + dot_v_n     * dt;
        next.x_g_n   = s[i].x_g_n   + dot_x_g_n   * dt;
        next.y_g_n   = s[i].y_g_n   + dot_y_g_n   * dt;

        next.Theta_r = s[i].Theta_r + dot_Theta_r * dt;
        next.v_r     = s[i].v_r     + dot_v_r     * dt;
        next.x_g_r   = s[i].x_g_r   + dot_x_g_r   * dt;
        next.y_g_r   = s[i].y_g_r   + dot_y_g_r   * dt;

        next.t = s[i].t + dt;

        s.push_back(next);
        ++i;
    }

    // --- Пост-обработка ---
    if (s.size() >= 2)
    {
        s[0].mass_r = s[1].mass_r;
        s[0].P_r = s[1].P_r;
    }

    std::cout << "Dist_to_target=" << s.back().r_rc
              << " Fuel_load=" << (m_fuel_raz + m_fuel_dva) <<  " beta=" << rocket.beta  << " mu=" << mu_0 << " k="<< MethodData.k << " kp=" << K_P << " km=" << K_m << " eta="<< eta_0 << std::endl;

    if (returnCode == "0")
    {
        std::cout << "Fuel_remain=" << m_fuel
                  << " m0_start=" << m_0_start << std::endl;
    }

    return returnCode;
}
