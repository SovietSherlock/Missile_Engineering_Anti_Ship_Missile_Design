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

double GetA(double M) {
    // Расчет значения вспомогательной функции A(M)
    double val = M * M - 1.0;
    double sign_val = static_cast<double>(sign(val));
    return  sign_val * std::sqrt(std::abs(val));
}

double GetSigmaM(double M) {
    // Расчет значения вспомогательной функции sigma(M)
    return 1.0 / (1.0 + std::exp(-M));
}

double GetCXBotCyl(double M) {
    // Расчет коэффициента донного сопротивления c_x_bot_cyl(M)


    // Промежуточные переменные A(M) и sigma(M)
    double A = GetA(M);
    double sig = GetSigmaM(M);

    // Система условий
    if (M <= 0.723672) {
        double inner = 1.410839 + 1.0432458 * A + 1.167756 * std::pow(A, 2.0) + 0.43818533 * std::pow(A, 3.0);
        return std::log(inner);
    }
    else if (M > 0.723672 && M <= 0.949985) {
        double power = -1.1994635 + 4.1244161 * A + 5.4425535 * std::pow(A, 2.0) + 2.6464556 * std::pow(A, 3.0);
        return std::exp(power);
    }
    else if (M > 0.949985 && M <= 1.045254) {
        double inner = 1508.25 - 6119.1935 * std::pow(M, -1.0) + 9304.258 * std::pow(M, -2.0) 
                       - 6277.8414 * std::pow(M, -3.0) + 1585.749 * std::pow(M, -4.0);
        return std::log(inner);
    }
    else if (M > 1.045254 && M <= 1.335822) {
        double power = -78.455781 + 249.63042 * M - 302.51885 * std::pow(M, 2.0) 
                       + 162.4847 * std::pow(M, 3.0) - 32.690523 * std::pow(M, 4.0);
        return std::exp(power);
    }
    else if (M > 1.335822 && M <= 3.74289) {
        double arg = 4.6685604 - 15.483104 * sig + 18.393009 * std::pow(sig, 2.0) - 7.5350128 * std::pow(sig, 3.0);
        return std::tan(arg);
    }
    else {
        double term1 = 1.43 / std::pow(M, 2.0);
        double term2 = (0.772 / std::pow(M, 2.0)) * std::pow(1.0 - 0.011 * std::pow(M, 2.0), 3.5);
        return term1 - term2;
    }
}

double GetCX (double cx0_pas, double M, double P, double S_a, double S_m) {
    // Расчет коэффициента лобового сопротивления 
    double cx0_bot_cyl = GetCXBotCyl(M);

    // 2. Предотвращаем деление на ноль, если площадь миделя передана некорректно
    if (S_m <= 0.0) return cx0_pas;

    // 3. Вычисляем компоненты по твоим встроенным функциям из LibConstFunc.h
    double sign_P = static_cast<double>(sign(P));     // Функция знака sgn(P)
    double heaviside_P = static_cast<double>(H(P));   // Функция Хевисайда H(P)

    // 4. Считаем по точной формуле: c_x = c_x_pas - c_x_0_дн * (S_a / S_m) * sgn(P) * H(P)
    double cx_total = cx0_pas - cx0_bot_cyl * (S_a / S_m) * sign_P * heaviside_P;

    // Физическая защита: полный коэффициент сопротивления не может стать отрицательным
    if (cx_total < 0.02) cx_total = 0.02;

    return cx_total;
};
