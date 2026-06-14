//Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.

#include <iostream>
#include <fstream>
#include <cmath>

#include "LibConstFunc.h"
#include "CalcTrajectory.h"
// #include "Zones.h"

using namespace std;

int main()
{
	setlocale(LC_CTYPE, "Russian");

	//Формируем структуры исходных данных для ПКР класса "Воздух-поверхность" с учетом установленных данных задания
	
	CalcParameters calcParams;	//Параметры расчёта
	calcParams.dt = 0.01;		//шаг по времени, с
	calcParams.r_por = 20.0;	//расстояние поражения, м

	InitConditions inits;			//Начальные условия

	// Параметры цели:
	inits.Theta_c0 = 0.0 / deg;		//угол направления вектора скорости цели, рад
	double v_c0_knots = 18.0; 		//скорость цели типа "Фрегат", м.узел
	inits.v_c0 = v_c0_knots * 0.514444;				//скорость цели, м/с
	inits.x_g_c0_min = 150.0*1e3;			//минимальная дальность полета ракеты до цели в условиях минимальной заданной высоты и скорости полета носителя, м
    inits.y_g_c0 = 0.0;			//координата цели, м

	// Параметры носителя:
	inits.Theta_n0 = 0.0 / deg;		//угол направления вектора скорости носителя, рад {если NaN, то равен Theta_r}	
	inits.v_n0 = 270.0;				//скорость носителя, м/с	
	inits.v_n0_min = 120.0;			//минимальная скорость носителя, м/с
	inits.x_g_pusk = 10.0*1e3;		//координата пуска по оси Y (начальная координата носителя и ракеты), м
	inits.y_g_pusk_min = 200.0;		//минимальная координата пуска по оси Y (минимальная начальная координата носителя и ракеты), м
	inits.y_g_pusk = 0.0;			//координата пуска по оси X (начальная координата носителя и ракеты), м

	// Параметры ракеты:
	inits.Theta_r0 = std::nan(""); //угол направления вектора скорости ракеты, рад {если NaN, то определяется методом наведения}	
    inits.v_r0 = 1000.0;				//скорость ракеты, м/с

	LimitConditions limits;			//Ограничения
    limits.n_ya_r_max = 40.0;       //максимальная перегрузка ракеты
	limits.t_r_min = 2.0;             //минимальное время полёта ракеты, с
	limits.t_r_max = 600.0;            //максимальное время полёта ракеты, с
	limits.R_r_min = nan("");         //минимальная дальность полёта ракеты, м {NaN: не проверяется}
	limits.R_r_max = nan("");         //максимальная дальность полёта ракеты, м {NaN: не проверяется}	
	limits.y_g_r_min = nan("");       //минимальная высота полёта ракеты, м {NaN: не проверяется}
	limits.y_g_r_max = nan("");       //максимальная высота полёта ракеты, м {NaN: не проверяется}
	limits.epsilon_nc_min = nan("") / deg;    //минимальный угол наклона линии визирования носитель-цель, рад {NaN: не проверяется}
	limits.epsilon_nc_max = nan("") / deg;    //максимальный угол наклона линии визирования носитель-цель, рад {NaN: не проверяется}

	C_Data CData;			//Данные по цели
	CData.n_xa_c_potr = {	//тангенциальная перегрузка цели (таблица время-перегрузка)
		{0.0, 0.0},
		{10.0, 0.0}
	};
	CData.n_ya_c_potr = {	//нормальная скоростная перегрузка цели  (таблица время-перегрузка)	
		{0.0, 1.0},
		{10.0, 1.0}
	};

	N_Data NData;			//Данные по носителю
	NData.n_xa_n_potr = {	//тангенциальная перегрузка носителя (таблица время-перегрузка)
		{0.0, 0.0},
		{10.0, 0.0}
	};
	NData.n_ya_n_potr = {	//нормальная скоростная перегрузка носителя  (таблица время-перегрузка)	
		{0.0, 1.0},
		{10.0, 1.0}
	};

	R_Data RData;
	RData.n_xa_r_potr = {	//тангенциальная перегрузка ракеты (таблица время-перегрузка)
		{0.0, 0.0},
		{10.0, 0.0}
	};

	GuidanceMethod_Data MethodData;		//Данные о методе наведения
	MethodData.Method = GuidanceMethod::BeamRiding_LineOfSight; //метод наведения (из перечисление доступных для выбора методов)
	MethodData.phi_upr = nan("") / deg; //угол упреждения, рад (для метода DeviatedPursuit, всегда положительный)
	MethodData.k = nan("");				//коэффициент пропорциональности (для метода ProportionalNav_const)
	MethodData.lambda = nan("");		//коэффициент пропорциональности (для метода ProportionalNav_var)
	MethodData.C = nan("");				//коэффициент в методе "С" (для метода BeamRiding_ConstAlignment)
	MethodData.m = nan("");				//коэффициент в методе спрямления (для метода BeamRiding_Alignment)

	//Объявляем структуру для резульата
	vector<TrajectoryParameters> result;

	string resultCode;
	std::ofstream C_File, N_File, R_File, Kinematic_File;
	C_File.open("Result_C.gra");
	N_File.open("Result_N.gra");
	R_File.open("Result_R.gra");
	Kinematic_File.open("Result_Kinematic.gra");
	
	//Вызываем функцию расчёта траектории
	resultCode = CalcTrajectory_Euler(calcParams, inits, limits, CData, NData, RData, MethodData, result);
	cout << "Результат: " << resultCode << endl;

	//Записываем результаты в файлы
    cout << "Запись файлов результатов..." << endl;

	C_File << "t\tTheta_c\tv_c\tx_g_c\ty_g_c\tn_xa_c\tn_ya_c";
	N_File << "t\tTheta_n\tv_n\tx_g_n\ty_g_n\tn_xa_n\tn_ya_n";
	R_File << "t\tTheta_r\tv_r\tx_g_r\ty_g_r\tn_xa_r\tn_ya_r";

	Kinematic_File << "t\tr_rc\tdot_r_rc\tddot_r_rc\tepsilon_rc\tdot_epsilon_rc\tddot_epsilon_rc" <<
		"\tr_nr\tdot_r_nr\tddot_r_nr\tepsilon_nr\tdot_epsilon_nr\tddot_epsilon_nr" <<
		"\tr_nc\tdot_r_nc\tddot_r_nc\tepsilon_nc\tdot_epsilon_nc\tddot_epsilon_nc" <<
		"\tDelta_r\tdot_Delta_r\tddot_Delta_r\tepsilon_l\tdot_epsilon_l\tddot_epsilon_l";

	for (size_t i = 0; i < result.size(); i++)
	{
		N_File << "\n" << result[i].t << "\t";
		N_File << result[i].Theta_n << "\t";
		N_File << result[i].v_n << "\t";
		N_File << result[i].x_g_n << "\t";
		N_File << result[i].y_g_n << "\t";
		N_File << result[i].n_xa_n << "\t";
		N_File << result[i].n_ya_n;

		R_File << "\n" << result[i].t << "\t";
		R_File << result[i].Theta_r << "\t";
		R_File << result[i].v_r << "\t";
		R_File << result[i].x_g_r << "\t";
		R_File << result[i].y_g_r << "\t";
		R_File << result[i].n_xa_r << "\t";
		R_File << result[i].n_ya_r;

		C_File << "\n" << result[i].t << "\t";
		C_File << result[i].Theta_c << "\t";
		C_File << result[i].v_c << "\t";
		C_File << result[i].x_g_c << "\t";
		C_File << result[i].y_g_c << "\t";
		C_File << result[i].n_xa_c << "\t";
		C_File << result[i].n_ya_c;

		Kinematic_File << "\n" << result[i].t << "\t";

		Kinematic_File << result[i].r_rc << "\t";
		Kinematic_File << result[i].dot_r_rc << "\t";
		Kinematic_File << result[i].ddot_r_rc << "\t";
		Kinematic_File << result[i].epsilon_rc << "\t";
		Kinematic_File << result[i].dot_epsilon_rc << "\t";
		Kinematic_File << result[i].ddot_epsilon_rc << "\t";

		Kinematic_File << result[i].r_nr << "\t";
		Kinematic_File << result[i].dot_r_nr << "\t";
		Kinematic_File << result[i].ddot_r_nr << "\t";
		Kinematic_File << result[i].epsilon_nr << "\t";
		Kinematic_File << result[i].dot_epsilon_nr << "\t";
		Kinematic_File << result[i].ddot_epsilon_nr << "\t";

		Kinematic_File << result[i].r_nc << "\t";
		Kinematic_File << result[i].dot_r_nc << "\t";
		Kinematic_File << result[i].ddot_r_nc << "\t";
		Kinematic_File << result[i].epsilon_nc << "\t";
		Kinematic_File << result[i].dot_epsilon_nc << "\t";
		Kinematic_File << result[i].ddot_epsilon_nc << "\t";

		Kinematic_File << result[i].Delta_r << "\t";
		Kinematic_File << result[i].dot_Delta_r << "\t";
		Kinematic_File << result[i].ddot_Delta_r << "\t";
		Kinematic_File << result[i].epsilon_l << "\t";
		Kinematic_File << result[i].dot_epsilon_l << "\t";
		Kinematic_File << result[i].ddot_epsilon_l;
	}

	N_File.close();
	C_File.close();
	R_File.close();
	Kinematic_File.close();

	cout << "Запись файлов завершена." << endl;	
		
	system("pause");	
}
