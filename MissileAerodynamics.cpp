#include "MissileAerodynamics.h"
#include "LibConstFunc.h" 
#include "Atmosphere.h"
#include <cmath>
#include <vector>


AeroCoefficients GetAerodynamics(double M) {
    AeroCoefficients aero;

    // Защита: число Маха не должно быть отрицательным
    if (M < 0.0) M = 0.0;

    // 6 опорных точек чисел Маха
    const double mach_points[] = { 0.4, 0.8, 1.1, 1.5, 3.0, 6.0 };
    
    // Значения коэффициентов
    const double cx0_points[]   = { 0.336, 0.348, 0.618, 0.571, 0.313, 0.170 };
    const double cya_points[]   = { 5.578, 10.768, 10.260, 9.487, 7.200, 5.741 };
    const double cyd_points[]   = { 5.548, 6.274, 7.233, 6.451, 3.432, 1.613 };
    const double bal_points[]   = { -0.098, -0.294, -0.210, -0.197, -0.204, -0.285 };

    // Упаковка данных в std::vector для твоей функции Linterp (размерность 6 строк на 2 столбца)
    std::vector<std::vector<double>> table_cx0(6, std::vector<double>(2));
    std::vector<std::vector<double>> table_cya(6, std::vector<double>(2));
    std::vector<std::vector<double>> table_cyd(6, std::vector<double>(2));
    std::vector<std::vector<double>> table_bal(6, std::vector<double>(2));

    for (int i = 0; i < 6; ++i) {
        table_cx0[i][0] = mach_points[i];
        table_cx0[i][1] = cx0_points[i];

        table_cya[i][0] = mach_points[i];
        table_cya[i][1] = cya_points[i];

        table_cyd[i][0] = mach_points[i];
        table_cyd[i][1] = cyd_points[i];

        table_bal[i][0] = mach_points[i];
        table_bal[i][1] = bal_points[i];
    }

    // Linterp автоматически найдет нужный отрезок Маха и посчитает точный коэффициент
    aero.cx0_pas         = Linterp(table_cx0, M, 1, 0);
    aero.cy_alpha        = Linterp(table_cya, M, 1, 0);
    aero.cy_delta        = Linterp(table_cyd, M, 1, 0);
    aero.alpha_delta_bal = Linterp(table_bal, M, 1, 0);

    return aero;
}

double GetEngineMass(){
    // Расчет массы топлива m_T
    double m_T;

    return m_T;
}

double GetKLambda(double d) {
    // Расчет коэффициента K_lambda(d)
    double denominator = std::pow(0.352105 * d, 3.044364);
    double exp_part = std::pow(0.9, -0.5 * d);
    double K_lambda = 4.388157e-5 * (1.0 / denominator) * exp_part;

    return K_lambda;
}

double GetEngineLambda(double d) {
    // Расчет длины двигательного отсека lambda(m_T, d)
    double K_lambda = GetKLambda(d);
    double m_T = GetEngineMass();
    return K_lambda * m_T;
}

double GetKBeta(double d){
    // Расчет коэффициента K_beta(d)
    double K_beta = 21.423397 * std::pow(1.605279 * d, 2.817663) + 0.079;
    return K_beta;
}

double GetBetaInf(double d){
    // Расчет коэффициента beta_infinity(d)
    double exp_part = std::pow(0.436906 * d, 0.107262);
    double beta_inf = 0.10013/exp_part + 1;
    return beta_inf;
}

double GetBeta(double d){
    // Расчет коэффициента beta(m_T, d)
    double K_beta = GetKBeta(d);
    double beta_inf = GetBetaInf(d);
    double m_T = GetEngineMass();
    double beta = K_beta * (1/m_T) + beta_inf;

}