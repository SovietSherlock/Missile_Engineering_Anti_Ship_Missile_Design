#pragma once

double n_ya_potr_PurePursuit(double dot_epsilon_rc, double v_r, double Theta_r);

double n_ya_potr_BeamRiding(double epsilon_l, double dot_epsilon_l, double ddot_epsilon_l,
                            double r_nr, double dot_r_nr, double Theta_r, double Theta_n,
                            double a_xa_r, double a_xa_n, double a_ya_n);

/// ПРОПОРЦИОНАЛЬНОЕ СБЛИЖЕНИЕ С ПОСТОЯННЫМ КОЭФФИЦИЕНТОМ
/// Формула (2.26): n_ya = (v_p/g) * k * dot_epsilon + cos(Theta_p)
double n_ya_potr_Proportional_const(double dot_epsilon_rc, double v_r, double k, double Theta_r);

/// ПРОПОРЦИОНАЛЬНОЕ СБЛИЖЕНИЕ С ПЕРЕМЕННЫМ КОЭФФИЦИЕНТОМ
/// Формула (2.41): n_ya = (lambda/g) * |dot_r| * dot_epsilon + cos(Theta_p)
double n_ya_potr_Proportional_var(double dot_epsilon_rc, double dot_r_rc, double lambda, double Theta_r);

double n_ya_potr_EVT(double v_r, double dot_epsilon_rc, double Theta_r, double K_g, double coeff);

double n_ya_potr_Gorka(double v_r, double Theta_r, double y_r, double r_rc,
                       double epsilon_rc, double h_prog_maneuver, double &dTheta_r_out);