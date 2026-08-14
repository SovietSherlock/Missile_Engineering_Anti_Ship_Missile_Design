#pragma once
// Определяем кинематические функции для расчёта параметров относительного движения двух точек

/// <summary>
/// Функция для вычисления расстояния между двумя точками
/// </summary>
/// <param name="x_g_r"> - координата x первой точки, м</param>
/// <param name="y_g_r"> - координата y первой точки, м</param>
/// <param name="x_g_c"> - координата x второй точки, м</param>
/// <param name="y_g_c"> - координата y второй точки, м</param>
/// <returns>Расстояние, м</returns>
double r(double x_g_r, double y_g_r, double x_g_c, double y_g_c);

/// <summary>
/// Функция для вычисления скорости изменения расстояния между двумя точками
/// </summary>
/// <param name="x_g_r"> - координата x первой точки, м</param>
/// <param name="y_g_r"> - координата y первой точки, м</param>
/// <param name="x_g_c"> - координата x второй точки, м</param>
/// <param name="y_g_c"> - координата y второй точки, м</param>
/// <param name="v_r"> - скорость первой точки, м/с</param>
/// <param name="v_c"> - скорость второй точки, м/с</param>
/// <param name="Theta_r"> - угол направления вектора скорости первой точки, рад</param>
/// <param name="Theta_c"> - угол направления вектора скорости второй точки, рад</param>
/// <returns>Скорость изменения расстояния, м/с</returns>
double dot_r(double x_g_r, double y_g_r, double x_g_c, double y_g_c, double v_r, double v_c, double Theta_r, double Theta_c);

/// <summary>
/// Функция для вычисления ускорения изменения расстояния между двумя точками
/// /// <summary>
/// <param name="x_g_r"> - координата x первой точки, м</param>
/// <param name="y_g_r"> - координата y первой точки, м</param>
/// <param name="x_g_c"> - координата x второй точки, м</param>
/// <param name="y_g_c"> - координата y второй точки, м</param>
/// <param name="v_r"> - скорость первой точки, м/с</param>
/// <param name="v_c"> - скорость второй точки, м/с</param>
/// <param name="Theta_r"> - угол направления вектора скорости первой точки, рад</param>
/// <param name="Theta_c"> - угол направления вектора скорости второй точки, рад</param>
/// <param name="a_xa_r"> - тангенциальное ускорение первой точки</param>
/// <param name="a_ya_r"> - нормальное скоростное ускорение первой точки</param>
/// <param name="a_xa_c"> - тангенциальное ускорение второй точки</param>
/// <param name="a_ya_c"> - нормальное скоростное ускорение второй точки</param>
/// <returns>Ускорение изменения расстояния, м/с^2</returns>
double ddot_r(double x_g_r, double y_g_r, double x_g_c, double y_g_c, double v_r, double v_c, double Theta_r, double Theta_c, double a_xa_r, double a_ya_r, double a_xa_c, double a_ya_c);

/// <summary>
/// Функция для вычисления угла наклона линии визирования между двумя точками
/// </summary>
/// <param name="x_g_r"> - координата x первой точки, м</param>
/// <param name="y_g_r"> - координата y первой точки, м</param>
/// <param name="x_g_c"> - координата x второй точки, м</param>
/// <param name="y_g_c"> - координата y второй точки, м</param>
/// <returns>Угол наклона линии визирования, рад</returns>
double epsilon(double x_g_r, double y_g_r, double x_g_c, double y_g_c);

/// <summary>
/// Функция для вычисления угловой скорости линии визирования между двумя точками
/// </summary>
/// <param name="x_g_r"> - координата x первой точки, м</param>
/// <param name="y_g_r"> - координата y первой точки, м</param>
/// <param name="x_g_c"> - координата x второй точки, м</param>
/// <param name="y_g_c"> - координата y второй точки, м</param>
/// <param name="v_r"> - скорость первой точки, м/с</param>
/// <param name="v_c"> - скорость второй точки, м/с</param>
/// <param name="Theta_r"> - угол направления вектора скорости первой точки, рад</param>
/// <param name="Theta_c"> - угол направления вектора скорости второй точки, рад</param>
/// <returns>Угловая скорость линии визирования, рад/с</returns>
double dot_epsilon(double x_g_r, double y_g_r, double x_g_c, double y_g_c, double v_r, double v_c, double Theta_r, double Theta_c);


/// <summary>
/// Функция для вычисления углового ускорения линии визирования между двумя точками
/// </summary>
/// <param name="x_g_r"> - координата x первой точки, м</param>
/// <param name="y_g_r"> - координата y первой точки, м</param>
/// <param name="x_g_c"> - координата x второй точки, м</param>
/// <param name="y_g_c"> - координата y второй точки, м</param>
/// <param name="v_r"> - скорость первой точки, м/с</param>
/// <param name="v_c"> - скорость второй точки, м/с</param>
/// <param name="Theta_r"> - угол направления вектора скорости первой точки, рад</param>
/// <param name="Theta_c"> - угол направления вектора скорости второй точки, рад</param>
/// <param name="a_xa_r"> - тангенциальное ускорение первой точки</param>
/// <param name="a_ya_r"> - нормальное скоростное ускорение первой точки</param>
/// <param name="a_xa_c"> - тангенциальное ускорение второй точки</param>
/// <param name="a_ya_c"> - нормальное скоростное ускорение второй точки</param>
/// <returns>Угловое ускорение линии визирования, рад/с^2</returns>
double ddot_epsilon(double x_g_r, double y_g_r, double x_g_c, double y_g_c, double v_r, double v_c, double Theta_r, double Theta_c, double a_xa_r, double a_ya_r, double a_xa_c, double a_ya_c);

/// <summary>
/// Функция для вычисления угловой скорости линии визирования между ракетой и ВИРТУАЛЬНОЙ целью
/// (используется для ЭВТ с виртуальной целью)
/// </summary>
/// <param name="x_g_r"> - координата x ракеты, м</param>
/// <param name="y_g_r"> - координата y ракеты, м</param>
/// <param name="x_g_c"> - координата x реальной цели, м</param>
/// <param name="y_g_c"> - координата y реальной цели, м</param>
/// <param name="delta_x_virtual"> - смещение виртуальной цели относительно реальной, м</param>
/// <param name="v_r"> - скорость ракеты, м/с</param>
/// <param name="v_c"> - скорость цели, м/с</param>
/// <param name="Theta_r"> - угол направления вектора скорости ракеты, рад</param>
/// <param name="Theta_c"> - угол направления вектора скорости цели, рад</param>
/// <returns>Угловая скорость линии визирования на виртуальную цель, рад/с</returns>
double dot_epsilon_virtual(
    double x_g_r, double y_g_r,
    double x_g_c, double y_g_c,
    double delta_x_virtual,
    double v_r, double v_c,
    double Theta_r, double Theta_c
);

/// <summary>
/// Функция для вычисления критического расстояния для начала пикирования
/// </summary>
/// <param name="v_r"> - скорость ракеты, м/с</param>
/// <param name="v_c"> - скорость цели, м/с</param>
/// <param name="H_gorka"> - высота горизонтального участка горки, м</param>
/// <param name="k_dive"> - коэффициент пропорциональной навигации для пикирования</param>
/// <param name="n_ya_max"> - максимальная располагаемая перегрузка</param>
/// <returns>Критическое расстояние для начала пикирования, м</returns>
double calc_r_critical(
    double v_r, double v_c,
    double H_gorka,
    double k_dive,
    double n_ya_max
);

/// <summary>
/// Функция для вычисления безопасной высоты полета ракеты по ЭВТ (границы радиогоризонта)
/// </summary>
/// <param name="r_rc"> - расстояние между ракетой и целью, м</param>
/// <param name="H_ant"> - высота антенны РЛС цели, м (по умолчанию 30)</param>
/// <param name="H_bez"> - высота безопасности, м (по умолчанию 0.5)</param>
/// <returns>Безопасная высота (граница радиогоризонта), м</returns>
double calc_radar_horizon_boundary(
    double r_rc,
    double H_ant = 30.0,
    double H_bez = 0.5
);