#pragma once

double n_ya_potr_PurePursuit(double dot_epsilon_rc, double v_r, double Theta_r);

double n_ya_potr_BeamRiding(double epsilon_l, double dot_epsilon_l, double ddot_epsilon_l, 
                            double r_nr, double dot_r_nr, double Theta_r, double Theta_n, 
                            double a_xa_r, double a_xa_n, double a_ya_n);

double n_ya_potr_Proportional(double k, double dot_epsilon_rc, double v_r, double Theta_r);

double n_ya_potr_EVT(double v_r, double dot_epsilon_rc, double Theta_r, double K_g, double coeff);

double n_ya_potr_EVT_virtual(double v_r, double v_c, double Theta_r,
                             double Theta_c, double dot_epsilon_virtual, double K_g, double k);

double n_ya_potr_March(double v_r, double Theta_r, double y_g_r, double H_march, double K_H_march, double K_v_march);

double n_ya_potr_GorkaClimb(double v_r, double Theta_r, double y_g_r,
                            double H_gorka, double K_H_gorka, double K_v_gorka);

double n_ya_potr_TwoPhase(
    double v_r, double Theta_r,
    double y_g_r,
    double r_rc,
    double epsilon_rc,
    double v_c, double Theta_c,
    double dot_epsilon_virtual,
    double dot_epsilon_rc,
    const GuidanceMethod_Data& mtd,
    int& phase,
    double r_crit_calc,
    double n_ya_max
);