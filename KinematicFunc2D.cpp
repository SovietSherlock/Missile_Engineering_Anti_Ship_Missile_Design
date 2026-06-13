#include "KinematicFunc2D.h"

#include <cmath>

#include "LibConstFunc.h"


/// Функция для вычисления расстояния между двумя точками
double r(double x_g_r, double y_g_r, double x_g_c, double y_g_c)
{
	return sqrt(sqr(x_g_c - x_g_r) + sqr(y_g_c - y_g_r));
}

/// Функция для вычисления скорости изменения расстояния между двумя точками
double dot_r(double x_g_r, double y_g_r, double x_g_c, double y_g_c, double v_r, double v_c, double Theta_r, double Theta_c)
{
	double r_rc = r(x_g_r, y_g_r, x_g_c, y_g_c);

	double Delta_x = x_g_c - x_g_r;
	double Delta_y = y_g_c - y_g_r;
	double Delta_dot_x = v_c * cos(Theta_c) - v_r * cos(Theta_r);
	double Delta_dot_y = v_c * sin(Theta_c) - v_r * sin(Theta_r);

	if (r_rc != 0.0)
		return (Delta_x * Delta_dot_x + Delta_y * Delta_dot_y) / r_rc;
	else
		return sqrt(sqr(Delta_dot_x) + sqr(Delta_dot_y));
}

/// Функция для вычисления ускорения изменения расстояния между двумя точками
double ddot_r(double x_g_r, double y_g_r, double x_g_c, double y_g_c, double v_r, double v_c, double Theta_r, double Theta_c, double a_xa_r, double a_ya_r, double a_xa_c, double a_ya_c)
{
	double r_rc = r(x_g_r, y_g_r, x_g_c, y_g_c);

	double dot_r_rc = dot_r(x_g_r, y_g_r, x_g_c, y_g_c, v_r, v_c, Theta_r, Theta_c);

	double Delta_x = x_g_c - x_g_r;
	double Delta_y = y_g_c - y_g_r;

	double Delta_dot_x = v_c * cos(Theta_c) - v_r * cos(Theta_r);
	double Delta_dot_y = v_c * sin(Theta_c) - v_r * sin(Theta_r);

	double Delta_ddot_x = a_xa_c * cos(Theta_c) - a_xa_r * cos(Theta_r) - a_ya_c * sin(Theta_c) + a_ya_r * sin(Theta_r);
	double Delta_ddot_y = a_xa_c * sin(Theta_c) - a_xa_r * sin(Theta_r) + a_ya_c * cos(Theta_c) - a_ya_r * cos(Theta_r);

	if (r_rc != 0.0)
		return (sqr(Delta_dot_x) + Delta_x * Delta_ddot_x + Delta_y * Delta_ddot_y + sqr(Delta_dot_y) - sqr(dot_r_rc)) / r_rc;
	else
		return sqrt(sqr(Delta_ddot_x) + sqr(Delta_ddot_y));
}

/// Функция для вычисления угла наклона линии визирования между двумя точками
double epsilon(double x_g_r, double y_g_r, double x_g_c, double y_g_c)
{
	if (x_g_c - x_g_r > 0.0)
		return atan((y_g_c - y_g_r) / (x_g_c - x_g_r));

	else if ((x_g_c - x_g_r < 0.0) && (y_g_c - y_g_r >= 0.0))
		return atan((y_g_c - y_g_r) / (x_g_c - x_g_r)) + pi;

	else if ((x_g_c - x_g_r < 0.0) && (y_g_c - y_g_r < 0.0))
		return atan((y_g_c - y_g_r) / (x_g_c - x_g_r)) - pi;

	else if ((x_g_c == x_g_r) && (y_g_c >= y_g_r))
		return pi / 2.0;

	else if ((x_g_c == x_g_r) && (y_g_c < y_g_r))
		return -pi / 2.0;
}

/// Функция для вычисления угловой скорости линии визирования между двумя точками
double dot_epsilon(double x_g_r, double y_g_r, double x_g_c, double y_g_c, double v_r, double v_c, double Theta_r, double Theta_c)
{
	double r_rc = r(x_g_r, y_g_r, x_g_c, y_g_c);

	double Delta_x = x_g_c - x_g_r;
	double Delta_y = y_g_c - y_g_r;

	double Delta_dot_x = v_c * cos(Theta_c) - v_r * cos(Theta_r);
	double Delta_dot_y = v_c * sin(Theta_c) - v_r * sin(Theta_r);

	return (Delta_x * Delta_dot_y - Delta_y * Delta_dot_x) / sqr(r_rc);
}

/// Функция для вычисления углового ускорения линии визирования между двумя точками
double ddot_epsilon(double x_g_r, double y_g_r, double x_g_c, double y_g_c, double v_r, double v_c, double Theta_r, double Theta_c, double a_xa_r, double a_ya_r, double a_xa_c, double a_ya_c)
{
	double r_rc = r(x_g_r, y_g_r, x_g_c, y_g_c);

	double Delta_x = x_g_c - x_g_r;
	double Delta_y = y_g_c - y_g_r;

	double Delta_dot_x = v_c * cos(Theta_c) - v_r * cos(Theta_r);
	double Delta_dot_y = v_c * sin(Theta_c) - v_r * sin(Theta_r);

	double Delta_ddot_x = a_xa_c * cos(Theta_c) - a_xa_r * cos(Theta_r) - a_ya_c * sin(Theta_c) + a_ya_r * sin(Theta_r);
	double Delta_ddot_y = a_xa_c * sin(Theta_c) - a_xa_r * sin(Theta_r) + a_ya_c * cos(Theta_c) - a_ya_r * cos(Theta_r);

	return (Delta_x * Delta_ddot_y - Delta_y * Delta_ddot_x) / sqr(r_rc)
		- 2.0 * (Delta_y * Delta_dot_y + Delta_x * Delta_dot_x) * (Delta_x * Delta_dot_y - Delta_y * Delta_dot_x) / pow(r_rc, 4.0);
}