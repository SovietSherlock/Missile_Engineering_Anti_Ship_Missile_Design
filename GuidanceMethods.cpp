#include "GuidanceMethods.h"

#include <cmath>

#include "LibConstFunc.h"

/// Функция для вычисления потребной перегрузки по методу чистой погони (самонаведение)
double n_ya_potr_PurePursuit(double dot_epsilon_rc, double v_r, double Theta_r)
{
    return v_r / g * dot_epsilon_rc + cos(Theta_r);
}

/// Функция для вычисления потребной перегрузки по методу наведения по лучу (теленаведение), для общего случая
double n_ya_potr_BeamRiding(double epsilon_l, double dot_epsilon_l, double ddot_epsilon_l, double r_nr, double dot_r_nr, double Theta_r, double Theta_n, double a_xa_r, double a_xa_n, double a_ya_n)
{
    return 1.0 / (g * cos(epsilon_l - Theta_r)) * (ddot_epsilon_l * r_nr + 2.0 * dot_epsilon_l * dot_r_nr +
                                                   a_xa_n * sin(Theta_n - epsilon_l) + a_ya_n * cos(epsilon_l - Theta_n) -
                                                   a_xa_r * sin(Theta_r - epsilon_l)) + cos(Theta_r);
}

double n_ya_potr_Proportional(double k, double dot_epsilon_rc, double v_r, double Theta_r)
{
    return v_r/ g * k * dot_epsilon_rc + cos(Theta_r);
}

double n_ya_potr_EVT(double v_r, double dot_epsilon_rc, double Theta_r, double K_g, double coeff)
{
    double result = (v_r / g) * dot_epsilon_rc * coeff + K_g * cos(Theta_r);
    return result;
}

double n_ya_potr_Gorka(double v_r, double Theta_r, double y_r, double r_rc,
                       double epsilon_rc, double h_prog_maneuver, double &dTheta_r_out)
{
    const double k_prog = 0.8;           // ↑ быстрее реакция, но без alpha_max
    const double h_cruise = 200.0;     // ↑ оптимум: ниже 300, выше 100
    const double h_peak = 300.0;       // подъём всего на +20 м
    const double theta_climb = 1.0;    // минимальный подъём для "горки"
    const double theta_dive_max = -35.0 * M_PI / 180.0; // не ограничиваем сильно

    dTheta_r_out = 0.0;
    double Theta_prog = 0.0;

    // === Фаза 1: Крейсерский полёт ===
    if (r_rc > h_prog_maneuver) {
        double h_err = h_cruise - y_r;
        Theta_prog = 0.002 * h_err;

        double theta_lim = 1.0 * M_PI / 180.0;
        if (Theta_prog >  theta_lim) Theta_prog =  theta_lim;
        if (Theta_prog < -theta_lim) Theta_prog = -theta_lim;

        dTheta_r_out = k_prog * (Theta_prog - Theta_r);
        double n_ya = (v_r * dTheta_r_out / g) + cos(Theta_r);
        if (n_ya > 1.1) n_ya = 1.1;
        if (n_ya < 0.9)  n_ya = 0.9;
        return n_ya;
    }

    // === Фаза 2: Параболический подъём (минимальный) ===
    if (r_rc > h_prog_maneuver * 0.25) {
        double t = (r_rc - h_prog_maneuver * 0.25) / (h_prog_maneuver * 0.75);
        double profile = std::sin(M_PI * t);  // плавная парабола
        Theta_prog = theta_climb * profile * M_PI / 180.0;

        // Жёсткий потолок
        if (y_r >= h_peak && Theta_prog > 0) Theta_prog = 0.0;

        dTheta_r_out = k_prog * (Theta_prog - Theta_r);
        double n_ya = (v_r * dTheta_r_out / g) + cos(Theta_r);
        if (n_ya > 1.3) n_ya = 1.3;
        if (n_ya < 0.8) n_ya = 0.8;
        return n_ya;
    }

    // === Фаза 3: Точное пикирование на цель ===
    Theta_prog = epsilon_rc;  // точно по линии визирования, без offset

    // Аварийное ограничение: не круче -35°
    if (Theta_prog < theta_dive_max) Theta_prog = theta_dive_max;

    // Если очень низко и не у цели — точно на цель
    if (y_r < 15.0 && r_rc > 100.0) Theta_prog = epsilon_rc;

    dTheta_r_out = k_prog * (Theta_prog - Theta_r);
    double n_ya = (v_r * dTheta_r_out / g) + cos(Theta_r);
    if (n_ya > 3.0) n_ya = 3.0;
    if (n_ya < -1.5) n_ya = -1.5;
    return n_ya;
}