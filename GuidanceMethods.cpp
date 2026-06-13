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
