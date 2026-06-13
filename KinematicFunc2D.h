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
