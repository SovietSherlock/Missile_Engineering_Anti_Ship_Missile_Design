#include "MissileAerodynamics.h"
#include "LibConstFunc.h" // Подключаем ради функции Linterp
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
