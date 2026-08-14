#pragma once
#include "DataStructs.h"
#include <vector>
#include <math.h>
#include <functional>

const double pi = 3.141592653589793; //число pi
const double deg = 57.295779513082320876798154814105; //отношение значения угла в градусах к значению в радианах
const double g = 9.80665; //ускорение свободного падения, м/с^2

#define sqr(x) ((x)*(x)) //функция квадрата

template <typename T>
int sign(T val) //функция сигнум (знак)
{
    return (T(0) < val) - (val < T(0));
}

template <typename T>
T H(T val) //функция Хевисайда
{
	if (val < 0.0)
		return 0.0;
	else if (val == 0.0)
		return 0.5;
	else
		return 1.0;
}

template <typename T>
T Linterp(const std::vector<std::vector<T>>& XY, T x, int Y_index = 1, int X_index = 0) //функция линейной интерполяции
{
    if (x <= XY[0][X_index])
    {
		return XY[0][Y_index];
	}
    else if(x >= XY[XY.size() - 1][X_index])
    {
		return XY[XY.size() - 1][Y_index];
	}
    else
    {
		size_t i = 0;
		while (i < XY.size() - 1)
		{
			if (x >= XY[i][X_index] && x < XY[i + 1][X_index])
			{
				return XY[i][Y_index] * ((XY[i + 1][X_index] - x) / (XY[i + 1][X_index] - XY[i][X_index])) + XY[i + 1][Y_index] * (1.0 - (XY[i + 1][X_index] - x) / (XY[i + 1][X_index] - XY[i][X_index]));
			}
			i = i + 1;
		}
	}
}

template <typename T>
T BiLinterp(const std::vector<std::vector<T>>& XYZ, T x, T y, int Z_index = 2) //функция билинейной интерполяции
{
	//ограничиваем по x
	T x_min = XYZ.front()[0];
	T x_max = XYZ.back()[0];
	if (x <= x_min)
	{
		size_t i = 0;
		T y_min = XYZ[i][1];
		T Z_pos;
		if (y <= y_min)
			Z_pos = XYZ[i][Z_index];
		else
		{
			while (XYZ[i][1] <= y && i < XYZ.size() - 1)
			{
				i++;
				if (XYZ[i - 1][1] > XYZ[i][1])
					break;
			}

			if (XYZ[i - 1][1] < XYZ[i][1])
				Z_pos = XYZ[i][Z_index] * ((XYZ[i - 1][1] - y) / (XYZ[i - 1][1] - XYZ[i][1])) + XYZ[i - 1][Z_index] * (1.0 - (XYZ[i - 1][1] - y) / (XYZ[i - 1][1] - XYZ[i][1]));
			else
				Z_pos = XYZ[i - 1][Z_index];
		}

		return Z_pos;
	}
	else if (x >= x_max)
	{
		size_t i = XYZ.size() - 1;
		T y_max = XYZ[i][1];
		T Z_pre;
		if (y >= y_max)
			Z_pre = XYZ[i][Z_index];
		else
		{
			while (XYZ[i][1] >= y && i > 0)
			{
				i--;
				if (XYZ[i][1] > XYZ[i + 1][1])
					break;
			}

			if (XYZ[i][1] < XYZ[i + 1][1])
				Z_pre = XYZ[i][Z_index] * ((XYZ[i + 1][1] - y) / (XYZ[i + 1][1] - XYZ[i][1])) + XYZ[i + 1][Z_index] * (1.0 - (XYZ[i + 1][1] - y) / (XYZ[i + 1][1] - XYZ[i][1]));
			else
				Z_pre = XYZ[i + 1][Z_index];
		}

		return Z_pre;
	}
	else
	{
		//ищем индекс начала области таблицы с первым x больше заданного
		size_t i_x = 0;
		while (XYZ[i_x][0] <= x && i_x < XYZ.size() - 1)
			i_x++;

		//ищем индекс y при меньшем х
		size_t i = i_x - 1;
		T y_max = XYZ[i][1];
		T Z_pre;
		if (y >= y_max)
			Z_pre = XYZ[i][Z_index];
		else
		{
			while (XYZ[i][1] >= y && i > 0)
			{
				i--;
				if (XYZ[i][1] > XYZ[i + 1][1])
					break;
			}

			if (XYZ[i][1] < XYZ[i + 1][1])
				Z_pre = XYZ[i][Z_index] * ((XYZ[i + 1][1] - y) / (XYZ[i + 1][1] - XYZ[i][1])) + XYZ[i + 1][Z_index] * (1.0 - (XYZ[i + 1][1] - y) / (XYZ[i + 1][1] - XYZ[i][1]));
			else
				Z_pre = XYZ[i + 1][Z_index];
		}

		//ищем индекс y при большем х
		i = i_x;
		T y_min = XYZ[i][1];
		T Z_pos;
		if (y <= y_min)
			Z_pos = XYZ[i][Z_index];
		else
		{
			while (XYZ[i][1] <= y && i < XYZ.size() - 1)
			{
				i++;
				if (XYZ[i - 1][1] > XYZ[i][1])
					break;
			}

			if (XYZ[i - 1][1] < XYZ[i][1])
				Z_pos = XYZ[i][Z_index] * ((XYZ[i - 1][1] - y) / (XYZ[i - 1][1] - XYZ[i][1])) + XYZ[i - 1][Z_index] * (1.0 - (XYZ[i - 1][1] - y) / (XYZ[i - 1][1] - XYZ[i][1]));
			else
				Z_pos = XYZ[i - 1][Z_index];
		}


		return Z_pre * ((XYZ[i_x][0] - x) / (XYZ[i_x][0] - XYZ[i_x - 1][0])) + Z_pos * (1.0 - (XYZ[i_x][0] - x) / (XYZ[i_x][0] - XYZ[i_x - 1][0]));
	}
}



inline double S_har(double d) {
    return pi * d * d / 4.0;
}

inline double S_a(double d) {
    return 0.8 * S_har(d);
}

// ============================================
// Линейная интерполяция
// ============================================
inline double linterp(double M, const double* Ms, const double* vals, int n)
{
    if (M <= Ms[0])     return vals[0];
    if (M >= Ms[n-1])   return vals[n-1];
    for (int i = 0; i < n - 1; ++i) {
        if (M >= Ms[i] && M <= Ms[i+1]) {
            double t = (M - Ms[i]) / (Ms[i+1] - Ms[i]);
            return vals[i] + t * (vals[i+1] - vals[i]);
        }
    }
    return vals[n-1];
}

// ============================================
// Таблица АД схемы
// ============================================
// inline double C_x0_pass(double M) {
//     static const double Ms[]   = {0.4, 0.8, 1.1, 1.5, 3.0, 6.0};
//     static const double vals[] = {0.511, 0.514, 0.832, 0.682, 0.4, 0.21};
//     return linterp(M, Ms, vals, 6);
// }

// inline double C_y_alpha(double M) {
//     static const double Ms[]   = {0.4, 0.8, 1.1, 1.5, 3.0, 6.0};
//     static const double vals[] = {18.427, 18.023, 21.8, 20.754, 16.864, 9.88};
//     return linterp(M, Ms, vals, 6);
// }

// inline double C_y_delta(double M) {
//     static const double Ms[]   = {0.4, 0.8, 1.1, 1.5, 3.0, 6.0};
//     static const double vals[] = {6.717, 7.44, 7.962, 6.722, 2.672, 0.947};
//     return linterp(M, Ms, vals, 6);
// }

// // (δ/α)_бал
// inline double delta_bal_over_alpha(double M) {
//     static const double Ms[]   = {0.4, 0.8, 1.1, 1.5, 3.0, 6.0};
//     static const double vals[] = {-0.404, -0.301, -0.348, -0.387, -0.78, -1.093};
//     return linterp(M, Ms, vals, 6);
// }

inline double C_x0_pass(double M)
{
    static const double Ms[]   = {0.4, 0.8, 1.1, 1.5, 3.0, 6.0};
    static const double vals[] = { 0.336, 0.348, 0.618, 0.571, 0.313, 0.170 };
    return linterp(M, Ms, vals, 6);
}

inline double C_y_alpha(double M)
{
    static const double Ms[]   = {0.4, 0.8, 1.1, 1.5, 3.0, 6.0};
    static const double vals[] = { 5.578, 10.768, 10.260, 9.487, 7.200, 5.741 };
    return linterp(M, Ms, vals, 6);
}

inline double C_y_delta(double M)
{
    static const double Ms[]   = {0.4, 0.8, 1.1, 1.5, 3.0, 6.0};
    static const double vals[] = { 5.548, 6.274, 7.233, 6.451, 3.432, 1.613 };
    return linterp(M, Ms, vals, 6);
}

inline double delta_bal_over_alpha(double M)
{
    static const double Ms[]   = {0.4, 0.8, 1.1, 1.5, 3.0, 6.0};
    static const double vals[] = { -0.098, -0.294, -0.210, -0.197, -0.204, -0.285 };
    return linterp(M, Ms, vals, 6);
}


// ============================================
// Донное сопротивление
// ============================================
inline double A_M(double M) {
    double mm1 = M * M - 1.0;
    return (mm1 >= 0.0 ? 1.0 : -1.0) * std::sqrt(std::fabs(mm1));
}

inline double sigma_M(double M) {
    return 1.0 / (1.0 + std::exp(-M));
}

inline double C_x0_don(double M) {
    if (M <= 0.723672) {
        double a = A_M(M);
        return std::log(1.410839 + 1.04324584*a + 1.167756*a*a + 0.43818533*a*a*a);
    }
    else if (M <= 0.949985) {
        double a = A_M(M);
        return std::exp(-1.1994635 + 4.1244161*a + 5.4425535*a*a + 2.6464556*a*a*a);
    }
    else if (M <= 1.045254) {
        double m1 = 1.0/M, m2 = m1*m1, m3 = m2*m1, m4 = m3*m1;
        return std::log(1508.25 - 6119.1935*m1 + 9304.258*m2 - 6277.8414*m3 + 1585.749*m4);
    }
    else if (M <= 1.335822) {
        double m2 = M*M, m3 = m2*M, m4 = m3*M;
        return std::exp(-78.455781 + 249.63042*M - 302.51885*m2 + 162.4847*m3 - 32.690523*m4);
    }
    else if (M <= 3.74289) {
        double s = sigma_M(M), s2 = s*s, s3 = s2*s;
        return std::tan(4.6685604 - 15.483104*s + 18.393009*s2 - 7.5350128*s3);
    }
    else {
        double m2 = M*M;
        double base = 1.0 - 0.011*m2;

        if (base <= 0.0)
        {
            return 1.43 / m2;
        }

        return (1.43 - 0.772 * std::pow(base, 3.5)) / m2;
    }
}

// ============================================
// Сборка C_x0 и C_y
// ============================================
inline double C_x0(double M, double d, double P) {
    double S_m = S_har(d);
    double sa  = S_a(d);
    double factor = (P > 0.0) ? 1.0 : ((P < 0.0) ? -1.0 : 0.0);   // sign(P)·H(P)
    return C_x0_pass(M) - C_x0_don(M) * (sa / S_m) * factor;
}

//  C_y = C_y^α·α + C_y^δ·δ
inline double C_y(double M, double alpha, double delta) {
    return C_y_alpha(M) * alpha + C_y_delta(M) * delta;
}

// ============================================
// Аэродинамические силы
// ============================================
inline double Y_a(double M, double alpha, double delta, double rho, double v, double d) {
    double S = S_har(d);
    return C_y(M, alpha, delta) * S * rho * v * v / 2.0;
}

// C_xa = C_x0 + C_y·sin(α)
// X_a  = C_xa · S · ρv²/2
inline double X_a(double M, double alpha, double delta, double rho, double v, double d, double P) {
    double S   = S_har(d);
    double cy  = C_y(M, alpha, delta);
    double cx0 = C_x0(M, d, P);
    double cxa = cx0 + cy * std::sin(alpha);
    return cxa * S * rho * v * v / 2.0;
}

// ============================================
// Поиск угла атаки α методом половинного деления
// ============================================
inline double find_alpha(double M, double rho, double v, double d,
                         double m, double P, double n_ya_req)
{
    if (v <= 0.0 || rho <= 0.0 || m <= 0.0 || d <= 0.0) {
        return 0.0;
    }

    const double tol      = 1e-6;
    const int    max_iter = 20;
    const double alpha_lim = 20.0 * M_PI / 180.0;

    double S      = S_har(d);
    double cy_a   = C_y_alpha(M);
    double cy_d   = C_y_delta(M);
    double dab    = delta_bal_over_alpha(M);
    double target = n_ya_req * m * 9.80665;

    // Если цель = 0 (прямолинейный полёт)
    if (std::abs(target) <= 0.00001f) {
        return 0.0;
    }

    auto f = [&](double alpha) -> double {
        double delta = dab * alpha;
        double cy    = cy_a * alpha + cy_d * delta;
        double Ya    = cy * S * rho * v * v / 2.0;
        return Ya + P * std::sin(alpha) - target;
    };

    // Интервал: от -alpha_lim до +alpha_lim
    double a = -alpha_lim;
    double b = alpha_lim;

    double fa = f(a);
    double fb = f(b);

    // Проверяем, есть ли корень в интервале
    if (fa * fb > 0.0) {
        // Корня нет — возвращаем границу с меньшим |f|
        double best = (std::fabs(fa) < std::fabs(fb)) ? a : b;
        return best;
    }

    // Метод половинного деления
    for (int i = 0; i < max_iter; ++i) {
        double c  = (a + b) / 2.0;
        double fc = f(c);
        if (std::fabs(fc) < tol || std::fabs(b - a) < tol)
            return c;
        if (fa * fc <= 0.0) {
            b = c; fb = fc;
        } else {
            a = c; fa = fc;
        }
    }
    return (a + b) / 2.0;
}

// Mu_0
inline double Mu_0(double m_start, double m_use_mass, double beta)
{
    return (m_start - m_use_mass) / (beta * m_start);
}

//beta
inline double calculate_beta(double d, double m_t)
{
    double K_d, beta_inf;

    if (d <= 0.165) {
        K_d = 3.502064 * std::pow(2.699929 * d, 2.2475);
        beta_inf = 0.010649 * (1.0 / std::pow(0.439973 * d, 0.962255)) + 1.0;
    } else {
        K_d = 21.423397 * std::pow(1.605279 * d, 2.817663) + 0.079;
        beta_inf = 0.10013 * (1.0 / std::pow(0.436906 * d, 0.107262)) + 1.0;
    }

    return K_d / m_t + beta_inf;
}