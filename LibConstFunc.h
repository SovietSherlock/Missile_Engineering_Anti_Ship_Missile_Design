#pragma once
#include <vector>
#include <cmath>
#include "DataStructs.h"

// ============================================================================
// Константы
// ============================================================================
const double pi = 3.141592653589793;
const double deg = 57.295779513082320876798154814105;
const double g = 9.80665;

#define sqr(x) ((x)*(x))

// ============================================================================
// Вспомогательные функции
// ============================================================================
template <typename T>
int sign(T val)
{
    return (T(0) < val) - (val < T(0));
}

template <typename T>
T H(T val)
{
    if (val < 0.0) return 0.0;
    if (val == 0.0) return 0.5;
    return 1.0;
}

// ============================================================================
// Линейная интерполяция
// ============================================================================
template <typename T>
T Linterp(const std::vector<std::vector<T>>& XY, T x, int Y_index = 1, int X_index = 0)
{
    if (x <= XY[0][X_index])
        return XY[0][Y_index];
    if (x >= XY[XY.size() - 1][X_index])
        return XY[XY.size() - 1][Y_index];

    for (size_t i = 0; i < XY.size() - 1; ++i)
    {
        if (x >= XY[i][X_index] && x < XY[i + 1][X_index])
        {
            double t = (x - XY[i][X_index]) / (XY[i + 1][X_index] - XY[i][X_index]);
            return XY[i][Y_index] * (1.0 - t) + XY[i + 1][Y_index] * t;
        }
    }
    return XY[XY.size() - 1][Y_index];
}

// ============================================================================
// Геометрия ракеты
// ============================================================================
inline double S_har(double d)
{
    return pi * d * d / 4.0;
}

inline double S_a(double d)
{
    return 0.8 * S_har(d);
}

// ============================================================================
// Линейная интерполяция для массивов
// ============================================================================
inline double linterp(double M, const double* Ms, const double* vals, int n)
{
    if (M <= Ms[0])     return vals[0];
    if (M >= Ms[n-1])   return vals[n-1];
    for (int i = 0; i < n - 1; ++i)
    {
        if (M >= Ms[i] && M <= Ms[i+1])
        {
            double t = (M - Ms[i]) / (Ms[i+1] - Ms[i]);
            return vals[i] + t * (vals[i+1] - vals[i]);
        }
    }
    return vals[n-1];
}

// ============================================================================
// Аэродинамические таблицы
// ============================================================================
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

// ============================================================================
// Донное сопротивление
// ============================================================================
inline double A_M(double M)
{
    double mm1 = M * M - 1.0;
    return (mm1 >= 0.0 ? 1.0 : -1.0) * std::sqrt(std::fabs(mm1));
}

inline double sigma_M(double M)
{
    return 1.0 / (1.0 + std::exp(-M));
}

inline double C_x0_don(double M)
{
    if (M <= 0.723672)
    {
        double a = A_M(M);
        return std::log(1.410839 + 1.04324584*a + 1.167756*a*a + 0.43818533*a*a*a);
    }
    else if (M <= 0.949985)
    {
        double a = A_M(M);
        return std::exp(-1.1994635 + 4.1244161*a + 5.4425535*a*a + 2.6464556*a*a*a);
    }
    else if (M <= 1.045254)
    {
        double m1 = 1.0/M, m2 = m1*m1, m3 = m2*m1, m4 = m3*m1;
        return std::log(1508.25 - 6119.1935*m1 + 9304.258*m2 - 6277.8414*m3 + 1585.749*m4);
    }
    else if (M <= 1.335822)
    {
        double m2 = M*M, m3 = m2*M, m4 = m3*M;
        return std::exp(-78.455781 + 249.63042*M - 302.51885*m2 + 162.4847*m3 - 32.690523*m4);
    }
    else if (M <= 3.74289)
    {
        double s = sigma_M(M), s2 = s*s, s3 = s2*s;
        return std::tan(4.6685604 - 15.483104*s + 18.393009*s2 - 7.5350128*s3);
    }
    else
    {
        double m2 = M*M;
        double base = 1.0 - 0.011*m2;
        if (base <= 0.0)
            return 1.43 / m2;
        return (1.43 - 0.772 * std::pow(base, 3.5)) / m2;
    }
}

// ============================================================================
// Сборка аэродинамических коэффициентов
// ============================================================================
inline double C_x0(double M, double d, double P)
{
    double S_m = S_har(d);
    double sa  = S_a(d);
    double factor = (P > 0.0) ? 1.0 : ((P < 0.0) ? -1.0 : 0.0);
    return C_x0_pass(M) - C_x0_don(M) * (sa / S_m) * factor;
}

inline double C_y(double M, double alpha, double delta)
{
    return C_y_alpha(M) * alpha + C_y_delta(M) * delta;
}

// ============================================================================
// Аэродинамические силы
// ============================================================================
inline double Y_a(double M, double alpha, double delta, double rho, double v, double d)
{
    double S = S_har(d);
    return C_y(M, alpha, delta) * S * rho * v * v / 2.0;
}

inline double X_a(double M, double alpha, double delta, double rho, double v, double d, double P)
{
    double S   = S_har(d);
    double cy  = C_y(M, alpha, delta);
    double cx0 = C_x0(M, d, P);
    double cxa = cx0 + cy * std::sin(alpha);
    return cxa * S * rho * v * v / 2.0;
}

// ============================================================================
// Поиск угла атаки α методом половинного деления
// ============================================================================
inline double find_alpha(double M, double rho, double v, double d,
                         double m, double P, double n_ya_req)
{
    if (v <= 0.0 || rho <= 0.0 || m <= 0.0 || d <= 0.0)
        return 0.0;

    const double tol      = 1e-6;
    const int    max_iter = 20;
    const double alpha_lim = 20.0 * pi / 180.0;

    double S      = S_har(d);
    double cy_a   = C_y_alpha(M);
    double cy_d   = C_y_delta(M);
    double dab    = delta_bal_over_alpha(M);
    double target = n_ya_req * m * 9.80665;

    if (std::fabs(target) < 1e-10)
        return 0.0;

    auto f = [&](double alpha) -> double {
        double delta = dab * alpha;
        double cy    = cy_a * alpha + cy_d * delta;
        double Ya    = cy * S * rho * v * v / 2.0;
        return Ya + P * std::sin(alpha) - target;
    };

    double a = -alpha_lim;
    double b = alpha_lim;
    double fa = f(a);
    double fb = f(b);

    if (fa * fb > 0.0)
    {
        double best = (std::fabs(fa) < std::fabs(fb)) ? a : b;
        return best;
    }

    for (int i = 0; i < max_iter; ++i)
    {
        double c  = (a + b) / 2.0;
        double fc = f(c);
        if (std::fabs(fc) < tol || std::fabs(b - a) < tol)
            return c;
        if (fa * fc <= 0.0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }
    return (a + b) / 2.0;
}

// ============================================================================
// Расчёт beta
// ============================================================================
inline double calculate_beta(double d, double m_t)
{
    double K_d, beta_inf;

    if (d <= 0.165)
    {
        K_d = 3.502064 * std::pow(2.699929 * d, 2.2475);
        beta_inf = 0.010649 * (1.0 / std::pow(0.439973 * d, 0.962255)) + 1.0;
    }
    else
    {
        K_d = 21.423397 * std::pow(1.605279 * d, 2.817663) + 0.079;
        beta_inf = 0.10013 * (1.0 / std::pow(0.436906 * d, 0.107262)) + 1.0;
    }

    return K_d / m_t + beta_inf;
}