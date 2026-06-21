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

// Функция вычисления массы топлива
double GetEngineMass();

// Функция для расчета коэффициента K_lambda в зависимости от диаметра d
double GetKLambda(double d);

// Функция для расчета удлинения/длины двигательной установки lambda(m_T, d)
double GetEngineLambda(double d);

// Функция вычисления коэффициента K_beta в зависимости от диаметра d
double GetKBeta(double d);

// Функция вычисления Beta_inf в зависимости от диаметра d
double GetBetaInf(double d);

// Функция вычисления Beta в зависимости от массы топлива m_T диаметра d
double GetBeta(double d);

#endif // MISSILEAERODYNAMICS_H
