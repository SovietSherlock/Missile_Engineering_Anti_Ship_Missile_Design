#ifndef MISSILEAERODYNAMICS_H
#define MISSILEAERODYNAMICS_H

// Структура для хранения мгновенных аэродинамических коэффициентов ракеты
struct AeroCoefficients {
    double cx0_pas;         // Коэффициент лобового сопротивления пассивный (при альфа=0, дельта=0)
    double cy_alpha;        // Производная коэффициента подъемной силы по углу атаки (1/рад)
    double cy_delta;        // Производная коэффициента подъемной силы по углу отклонения рулей (1/рад)
    double alpha_delta_bal; // Балансировочное соотношение (alpha / delta)_бал
};

// Функция, вычисляющая коэффициенты для заданного числа Маха (M)
AeroCoefficients GetAerodynamics(double M);

#endif // MISSILEAERODYNAMICS_H
