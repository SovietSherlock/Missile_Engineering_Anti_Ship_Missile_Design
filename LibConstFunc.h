#pragma once
#include <vector>

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
