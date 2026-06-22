#include "GuidanceMethods.h"
#include <cmath>
#include "LibConstFunc.h"

double n_ya_potr_PurePursuit(double dot_epsilon_rc, double v_r, double Theta_r)
{
    return v_r / g * dot_epsilon_rc + cos(Theta_r);
}

double n_ya_potr_BeamRiding(double epsilon_l, double dot_epsilon_l, double ddot_epsilon_l,
                            double r_nr, double dot_r_nr, double Theta_r, double Theta_n,
                            double a_xa_r, double a_xa_n, double a_ya_n)
{
    return 1.0 / (g * cos(epsilon_l - Theta_r)) *
           (ddot_epsilon_l * r_nr + 2.0 * dot_epsilon_l * dot_r_nr +
            a_xa_n * sin(Theta_n - epsilon_l) + a_ya_n * cos(epsilon_l - Theta_n) -
            a_xa_r * sin(Theta_r - epsilon_l)) + cos(Theta_r);
}

double n_ya_potr_Proportional_const(double dot_epsilon_rc, double v_r, double k, double Theta_r)
{
    return (v_r / g) * k * dot_epsilon_rc + cos(Theta_r);
}

double n_ya_potr_Proportional_var(double dot_epsilon_rc, double dot_r_rc, double lambda, double Theta_r)
{
    return (lambda / g) * fabs(dot_r_rc) * dot_epsilon_rc + cos(Theta_r);
}

double n_ya_potr_EVT(double v_r, double dot_epsilon_rc, double Theta_r, double K_g, double coeff)
{
    return (v_r / g) * dot_epsilon_rc * coeff + K_g * cos(Theta_r);
}

double n_ya_potr_Gorka(double v_r, double Theta_r, double y_r, double r_rc,
                       double epsilon_rc, double h_prog_maneuver, double &dTheta_r_out)
{
    const double k_prog = 0.8;
    const double h_cruise = 200.0;
    const double h_peak = 300.0;
    const double theta_climb = 1.0 * pi / 180.0;
    const double theta_dive_max = -35.0 * pi / 180.0;

    dTheta_r_out = 0.0;
    double Theta_prog = 0.0;

    // Фаза 1: Крейсерский полёт
    if (r_rc > h_prog_maneuver)
    {
        double h_err = h_cruise - y_r;
        Theta_prog = 0.002 * h_err;

        double theta_lim = 1.0 * pi / 180.0;
        if (Theta_prog >  theta_lim) Theta_prog =  theta_lim;
        if (Theta_prog < -theta_lim) Theta_prog = -theta_lim;

        dTheta_r_out = k_prog * (Theta_prog - Theta_r);
        double n_ya = (v_r * dTheta_r_out / g) + cos(Theta_r);
        if (n_ya > 1.1) n_ya = 1.1;
        if (n_ya < 0.9) n_ya = 0.9;
        return n_ya;
    }

    // Фаза 2: Параболический подъём
    if (r_rc > h_prog_maneuver * 0.25)
    {
        double t = (r_rc - h_prog_maneuver * 0.25) / (h_prog_maneuver * 0.75);
        double profile = std::sin(pi * t);
        Theta_prog = theta_climb * profile;

        if (y_r >= h_peak && Theta_prog > 0) Theta_prog = 0.0;

        dTheta_r_out = k_prog * (Theta_prog - Theta_r);
        double n_ya = (v_r * dTheta_r_out / g) + cos(Theta_r);
        if (n_ya > 1.3) n_ya = 1.3;
        if (n_ya < 0.8) n_ya = 0.8;
        return n_ya;
    }

    // Фаза 3: Пикирование
    Theta_prog = epsilon_rc;
    if (Theta_prog < theta_dive_max) Theta_prog = theta_dive_max;
    if (y_r < 15.0 && r_rc > 100.0) Theta_prog = epsilon_rc;

    dTheta_r_out = k_prog * (Theta_prog - Theta_r);
    double n_ya = (v_r * dTheta_r_out / g) + cos(Theta_r);
    if (n_ya > 3.0) n_ya = 3.0;
    if (n_ya < -1.5) n_ya = -1.5;
    return n_ya;
}