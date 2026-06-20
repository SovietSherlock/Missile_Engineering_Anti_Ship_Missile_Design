#include "Atmosphere.h"
#include <cmath>

AirData Atmosphere_GOST_4401_81(double h){
    AirData state;

    // Основные физические константы:
    const double g_0 = 9.80665;
    const double r_earth = 6356767.0;
    const double R = 8314.32 / 28.96442;
    const double beta_S = 1.458e-6;
    const double S = 110.4;
    const double k = 1.4;

    if (h < 0.0) h = 0.0;

    // Расчет ускорения свободного падения g(h)
    state.g = g_0 * std::pow((r_earth / (r_earth + h)), 2);

    // Перевод геометрической высоты (h) в геопотенциальную (H)
    double H = r_earth * h / (r_earth + h);

    // Массивы опорных точек ГОСТ 4401-81 из твоего скрипта
    // Индексы слоев: 0: 0-11км, 1: 11-20км, 2: 20-32км, 3: 32-47км, 4: 47-51км, 5: 51-71км, 6: 71-85км
    const double data_H[] = { 0.0, 11000.0, 20000.0, 32000.0, 47000.0, 51000.0, 71000.0, 85000.0 };
    const double data_T[] = { 288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65, 186.65 };
    const double data_B[] = { -6.5e-3, 0.0, 1.0e-3, 2.8e-3, 0.0, -2.8e-3, -2.0e-3 };
    const double data_p[] = { 101325.0, 22632.06, 5474.889, 868.0187, 110.9063, 66.93887, 3.956420, 0.3734 }; 
    // Примечание: p_data рассчитаны заранее по формулам твоего __init__ для ускорения работы C++

    // Определение текущего высотного слоя
    int i = 0;
    if (H < data_H[0]) {
        state.T = data_T[0];
        state.p = data_p[0];
    } 
    else if (H >= data_H[7]) {
        state.T = data_T[7];
        state.p = data_p[7];
    } 
    else {
        // Ищем, в какой интервал попала высота H
        for (int j = 0; j < 7; ++j) {
            if (H >= data_H[j]) {
                i = j;
            }
        }

        // Расчет температуры T(H) через линейный градиент (аналог interp1d)
        state.T = data_T[i] + data_B[i] * (H - data_H[i]);

        // Расчет давления p(H) в зависимости от градиента B
        if (data_B[i] != 0.0) {
            state.p = data_p[i] * std::pow((state.T / data_T[i]), (-g_0 / (data_B[i] * R)));
        } else {
            state.p = data_p[i] * std::exp(-g_0 * (H - data_H[i]) / (R * data_T[i]));
        }
    }

    // Расчет плотности rho(h)
    state.rho = state.p / (R * state.T);

    // Расчет скорости звука a(h)
    state.a = std::sqrt(k * R * state.T);

    // Расчет динамической вязкости mu(T) по Сазерленду
    state.mu = beta_S * std::pow(state.T, 1.5) / (state.T + S);

    // Расчет теплопроводности lam(T)
    state.lam = 2.648151e-3 * std::pow(state.T, 1.5) / (state.T + 245.4 * std::pow(10.0, (-12.0 / state.T)));

    return state;
}