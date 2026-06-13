#pragma once
// Определяем функции для вычисления потребной перегрузки по различным методам наведения

/// Функция для вычисления потребной перегрузки по методу чистой погони (самонаведение)
double n_ya_potr_PurePursuit(double dot_epsilon_rc, double v_r, double Theta_r);

/// Функция для вычисления потребной перегрузки по методу совмещения (теленаведение), для общего случая
double n_ya_potr_BeamRiding(double epsilon_l, double dot_epsilon_l, double ddot_epsilon_l, double r_nr, double dot_r_nr, double Theta_r, double Theta_n, double a_xa_r, double a_xa_n, double a_ya_n);