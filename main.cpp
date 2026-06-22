#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <cfloat>
#include <chrono>
#include <algorithm>
#include <functional>
#include <random>
#include <map>
#include <tuple>
#include <sys/stat.h>
#include <direct.h>
#include "DataStructs.h"
#include "CalcTrajectory.h"
#include "LibConstFunc.h"

// ============================================================================
// Константы
// ============================================================================
const double PI = 3.14159265358979323846;

// ============================================================================
// Создание папки
// ============================================================================
void createDirectory(const std::string& path)
{
#ifdef _WIN32
    _mkdir(path.c_str());
#else
    mkdir(path.c_str(), 0755);
#endif
}

// ============================================================================
// Вспомогательная функция для вывода в файл
// ============================================================================
void saveTrajectoryToFile(const std::vector<TrajectoryParameters>& s,
                          const std::string& filename,
                          const std::string& header,
                          const std::vector<std::string>& fields)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    file << header << "\n";
    for (const auto& p : s)
    {
        file << std::setprecision(10) << std::scientific;
        file << p.t;

        for (const auto& field : fields)
        {
            if (field == "Theta_c") file << "\t" << p.Theta_c;
            else if (field == "v_c") file << "\t" << p.v_c;
            else if (field == "x_g_c") file << "\t" << p.x_g_c;
            else if (field == "y_g_c") file << "\t" << p.y_g_c;
            else if (field == "n_xa_c") file << "\t" << p.n_xa_c;
            else if (field == "n_ya_c") file << "\t" << p.n_ya_c;
            else if (field == "Theta_n") file << "\t" << p.Theta_n;
            else if (field == "v_n") file << "\t" << p.v_n;
            else if (field == "x_g_n") file << "\t" << p.x_g_n;
            else if (field == "y_g_n") file << "\t" << p.y_g_n;
            else if (field == "n_xa_n") file << "\t" << p.n_xa_n;
            else if (field == "n_ya_n") file << "\t" << p.n_ya_n;
            else if (field == "Theta_r") file << "\t" << p.Theta_r;
            else if (field == "v_r") file << "\t" << p.v_r;
            else if (field == "x_g_r") file << "\t" << p.x_g_r;
            else if (field == "y_g_r") file << "\t" << p.y_g_r;
            else if (field == "n_xa_r") file << "\t" << p.n_xa_r;
            else if (field == "n_ya_r") file << "\t" << p.n_ya_r;
            else if (field == "r_rc") file << "\t" << p.r_rc;
            else if (field == "dot_r_rc") file << "\t" << p.dot_r_rc;
            else if (field == "ddot_r_rc") file << "\t" << p.ddot_r_rc;
            else if (field == "epsilon_rc") file << "\t" << p.epsilon_rc;
            else if (field == "dot_epsilon_rc") file << "\t" << p.dot_epsilon_rc;
            else if (field == "ddot_epsilon_rc") file << "\t" << p.ddot_epsilon_rc;
            else if (field == "mass_r") file << "\t" << p.mass_r;
            else if (field == "P_r") file << "\t" << p.P_r;
            else if (field == "alpha_r") file << "\t" << p.alpha_r;
            else if (field == "delta_r") file << "\t" << p.delta_r;
            else if (field == "mode1") file << "\t" << (p.mode1 ? 1 : 0);
        }
        file << "\n";
    }
    file.close();
    std::cout << "  Saved: " << filename << std::endl;
}

// ============================================================================
// Функция для сохранения данных графиков
// ============================================================================
void savePlotData(const std::vector<TrajectoryParameters>& s,
                  const std::string& filename,
                  const std::string& title,
                  const std::vector<std::pair<std::string, std::function<double(const TrajectoryParameters&)>>>& columns)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    file << "# " << title << std::endl;
    for (const auto& col : columns)
    {
        file << "# " << col.first;
    }
    file << std::endl;
    file << std::setprecision(10) << std::scientific;

    for (const auto& p : s)
    {
        for (size_t i = 0; i < columns.size(); ++i)
        {
            if (i > 0) file << "\t";
            file << columns[i].second(p);
        }
        file << "\n";
    }
    file.close();
    std::cout << "  Saved plot: " << filename << std::endl;
}

// ============================================================================
// Функции для графиков
// ============================================================================
void saveTrajectoryPlot(const std::vector<TrajectoryParameters>& s,
                        const std::string& filename,
                        const std::string& title)
{
    std::vector<std::pair<std::string, std::function<double(const TrajectoryParameters&)>>> cols = {
        {"x_g", [](const TrajectoryParameters& p) { return p.x_g_r; }},
        {"y_g", [](const TrajectoryParameters& p) { return p.y_g_r; }}
    };
    savePlotData(s, filename, title, cols);
}

void saveVelocityPlot(const std::vector<TrajectoryParameters>& s,
                      const std::string& filename,
                      const std::string& title)
{
    std::vector<std::pair<std::string, std::function<double(const TrajectoryParameters&)>>> cols = {
        {"t", [](const TrajectoryParameters& p) { return p.t; }},
        {"v", [](const TrajectoryParameters& p) { return p.v_r; }}
    };
    savePlotData(s, filename, title, cols);
}

void saveCoordPlots(const std::vector<TrajectoryParameters>& s,
                    const std::string& prefix,
                    const std::string& title)
{
    std::vector<std::pair<std::string, std::function<double(const TrajectoryParameters&)>>> cols_x = {
        {"t", [](const TrajectoryParameters& p) { return p.t; }},
        {"x_g", [](const TrajectoryParameters& p) { return p.x_g_r; }}
    };
    savePlotData(s, prefix + "_x_t.gra", title + " - x_g(t)", cols_x);

    std::vector<std::pair<std::string, std::function<double(const TrajectoryParameters&)>>> cols_y = {
        {"t", [](const TrajectoryParameters& p) { return p.t; }},
        {"y_g", [](const TrajectoryParameters& p) { return p.y_g_r; }}
    };
    savePlotData(s, prefix + "_y_t.gra", title + " - y_g(t)", cols_y);
}

void saveThetaPlot(const std::vector<TrajectoryParameters>& s,
                   const std::string& filename,
                   const std::string& title)
{
    std::vector<std::pair<std::string, std::function<double(const TrajectoryParameters&)>>> cols = {
        {"t", [](const TrajectoryParameters& p) { return p.t; }},
        {"Theta_rad", [](const TrajectoryParameters& p) { return p.Theta_r; }},
        {"Theta_deg", [](const TrajectoryParameters& p) { return p.Theta_r * 180.0 / PI; }}
    };
    savePlotData(s, filename, title, cols);
}

void saveOverloadPlot(const std::vector<TrajectoryParameters>& s,
                      const std::string& filename,
                      const std::string& title)
{
    std::vector<std::pair<std::string, std::function<double(const TrajectoryParameters&)>>> cols = {
        {"t", [](const TrajectoryParameters& p) { return p.t; }},
        {"n_ya", [](const TrajectoryParameters& p) { return p.n_ya_r; }}
    };
    savePlotData(s, filename, title, cols);
}

void saveAltitudePlot(const std::vector<TrajectoryParameters>& s,
                      const std::string& filename,
                      double y_opt,
                      const std::string& title)
{
    std::vector<std::pair<std::string, std::function<double(const TrajectoryParameters&)>>> cols = {
        {"t", [](const TrajectoryParameters& p) { return p.t; }},
        {"y_g", [](const TrajectoryParameters& p) { return p.y_g_r; }},
        {"y_opt", [y_opt](const TrajectoryParameters&) { return y_opt; }}
    };
    savePlotData(s, filename, title, cols);
}

// ============================================================================
// Gnuplot-скрипт
// ============================================================================
void saveGnuplotScript(const std::string& filename,
                       const std::string& title,
                       const std::string& prefix)
{
    std::ofstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    file << "# Gnuplot script for " << title << std::endl;
    file << "# Usage: gnuplot -persist " << filename << std::endl;
    file << std::endl;

    file << "set terminal wxt size 800,600 enhanced font 'Arial,12'" << std::endl;
    file << "set title '" << title << " - Траектория полёта'" << std::endl;
    file << "set xlabel 'x_g, м'" << std::endl;
    file << "set ylabel 'y_g, м'" << std::endl;
    file << "set grid" << std::endl;
    file << "plot '" << prefix << "_trajectory.gra' using 1:2 with lines lw 2 lt rgb 'red' title 'Ракета'" << std::endl;
    file << "pause -1 'Нажмите Enter для следующего графика...'" << std::endl;
    file << std::endl;

    file << "set title '" << title << " - Скорость ракеты'" << std::endl;
    file << "set xlabel 't, с'" << std::endl;
    file << "set ylabel 'v, м/с'" << std::endl;
    file << "set grid" << std::endl;
    file << "plot '" << prefix << "_velocity.gra' using 1:2 with lines lw 2 lt rgb 'blue' title 'v(t)'" << std::endl;
    file << "pause -1 'Нажмите Enter для следующего графика...'" << std::endl;
    file << std::endl;

    file << "set title '" << title << " - Координаты'" << std::endl;
    file << "set xlabel 't, с'" << std::endl;
    file << "set ylabel 'Координата, м'" << std::endl;
    file << "set grid" << std::endl;
    file << "plot '" << prefix << "_x_t.gra' using 1:2 with lines lw 2 lt rgb 'green' title 'x_g(t)', \\" << std::endl;
    file << "     '" << prefix << "_y_t.gra' using 1:2 with lines lw 2 lt rgb 'red' title 'y_g(t)'" << std::endl;
    file << "pause -1 'Нажмите Enter для следующего графика...'" << std::endl;
    file << std::endl;

    file << "set title '" << title << " - Угол траектории'" << std::endl;
    file << "set xlabel 't, с'" << std::endl;
    file << "set ylabel 'Theta, град'" << std::endl;
    file << "set grid" << std::endl;
    file << "plot '" << prefix << "_theta.gra' using 1:3 with lines lw 2 lt rgb 'purple' title 'Theta(t)'" << std::endl;
    file << "pause -1 'Нажмите Enter для следующего графика...'" << std::endl;
    file << std::endl;

    file << "set title '" << title << " - Нормальная перегрузка'" << std::endl;
    file << "set xlabel 't, с'" << std::endl;
    file << "set ylabel 'n_ya'" << std::endl;
    file << "set grid" << std::endl;
    file << "plot '" << prefix << "_n_ya.gra' using 1:2 with lines lw 2 lt rgb 'orange' title 'n_ya(t)'" << std::endl;
    file << "pause -1 'Нажмите Enter для завершения...'" << std::endl;

    file.close();
    std::cout << "  Saved Gnuplot script: " << filename << std::endl;
}

// ============================================================================
// Функция для сохранения ВСЕХ графиков в папку
// ============================================================================
void saveAllPlots(const std::vector<TrajectoryParameters>& s,
                  const std::string& prefix,
                  const std::string& label,
                  double y_opt = 10000.0)
{
    if (s.empty())
    {
        std::cout << "  Нет данных для построения графиков" << std::endl;
        return;
    }

    std::string plotDir = "plots";
    createDirectory(plotDir);

    std::cout << "\n--- ГЕНЕРАЦИЯ ГРАФИКОВ: " << label << " ---" << std::endl;
    std::cout << "  Папка: " << plotDir << "/" << std::endl;

    auto getPath = [&](const std::string& name) -> std::string {
        return plotDir + "/" + prefix + "_" + name;
    };

    saveTrajectoryPlot(s, getPath("trajectory.gra"), label + ": y_g(x_g)");
    saveVelocityPlot(s, getPath("velocity.gra"), label + ": v(t)");
    saveCoordPlots(s, getPath(""), label);
    saveThetaPlot(s, getPath("theta.gra"), label + ": Theta(t)");
    saveOverloadPlot(s, getPath("n_ya.gra"), label + ": n_ya(t)");
    saveAltitudePlot(s, getPath("altitude.gra"), y_opt, label + ": Высота полёта");
    saveGnuplotScript(getPath("plot.gnu"), label, prefix);

    std::cout << "  Все графики сохранены в папку: " << plotDir << "/" << std::endl;
}

// ============================================================================
// Функция для сохранения результатов
// ============================================================================
void saveAllResults(const std::vector<TrajectoryParameters>& s,
                    const std::string& prefix)
{
    if (s.empty())
    {
        std::cout << "  Нет данных для сохранения" << std::endl;
        return;
    }

    std::cout << "\nСохранение результатов с префиксом: " << prefix << std::endl;

    std::vector<std::string> cFields = {"Theta_c", "v_c", "x_g_c", "y_g_c", "n_xa_c", "n_ya_c"};
    saveTrajectoryToFile(s, prefix + "_C.gra",
                         "t\tTheta_c\tv_c\tx_g_c\ty_g_c\tn_xa_c\tn_ya_c",
                         cFields);

    std::vector<std::string> nFields = {"Theta_n", "v_n", "x_g_n", "y_g_n", "n_xa_n", "n_ya_n"};
    saveTrajectoryToFile(s, prefix + "_N.gra",
                         "t\tTheta_n\tv_n\tx_g_n\ty_g_n\tn_xa_n\tn_ya_n",
                         nFields);

    std::vector<std::string> rFields = {"Theta_r", "v_r", "x_g_r", "y_g_r", "n_xa_r", "n_ya_r",
                                        "mass_r", "P_r", "alpha_r", "delta_r", "mode1"};
    saveTrajectoryToFile(s, prefix + "_R.gra",
                         "t\tTheta_r\tv_r\tx_g_r\ty_g_r\tn_xa_r\tn_ya_r\tmass_r\tP_r\talpha_r\tdelta_r\tmode1",
                         rFields);

    std::vector<std::string> kFields = {"r_rc", "dot_r_rc", "ddot_r_rc", "epsilon_rc",
                                        "dot_epsilon_rc", "ddot_epsilon_rc"};
    saveTrajectoryToFile(s, prefix + "_Kinematic.gra",
                         "t\tr_rc\tdot_r_rc\tddot_r_rc\tepsilon_rc\tdot_epsilon_rc\tddot_epsilon_rc",
                         kFields);

    std::cout << "  Сохранено " << s.size() << " точек" << std::endl;
}

// ============================================================================
// ФУНКЦИЯ ДЛЯ ВЫПОЛНЕНИЯ РАСЧЁТА С ЗАДАННЫМИ ПАРАМЕТРАМИ
// ============================================================================
struct TrajectoryResult
{
    double range_km;
    double flight_time;
    std::string result_code;
    std::vector<TrajectoryParameters> trajectory;
    bool is_success;
};

TrajectoryResult runTrajectoryCalculation(
    const CalcParameters& calcParams,
    const InitConditions& inits,
    const LimitConditions& limits,
    const C_Data& CData,
    const N_Data& NData,
    const R_Data& RData,
    const RocketParams& rocket,
    double mu,
    double eta,
    double K_m,
    double K_P,
    double K_g,
    double m_norm,
    double switch_time)
{
    TrajectoryResult result;
    result.is_success = false;
    result.range_km = 0.0;
    result.flight_time = 0.0;

    GuidanceMethod_Data methodData;
    methodData.Method = GuidanceMethod::EVT;
    methodData.k = 3.0;
    methodData.K_g = K_g;

    std::vector<TrajectoryParameters> s;
    std::string resCode = CalcTrajectory_DZ(
        calcParams, inits, limits,
        CData, NData, RData,
        methodData,
        rocket,
        mu, eta, K_m, K_P, K_g,
        s,
        m_norm,
        switch_time
    );

    result.trajectory = s;
    result.result_code = resCode;

    if (resCode == "0" && !s.empty())
    {
        result.is_success = true;
        result.flight_time = s.back().t;
        double x_start = s.front().x_g_r;
        double x_end = s.back().x_g_r;
        result.range_km = (x_end - x_start) / 1000.0;
    }

    return result;
}

// ============================================================================
// ФУНКЦИЯ ДЛЯ СКАНИРОВАНИЯ ПАРАМЕТРОВ И ПОИСКА МАКСИМУМА ДЛЯ ВСЕХ МЕТОДОВ
// ============================================================================
void scanAndOptimizeAllMethods(
    const CalcParameters& calcParams,
    const InitConditions& inits,
    const LimitConditions& limits,
    const C_Data& CData,
    const N_Data& NData,
    const R_Data& RData,
    const RocketParams& rocket,
    double K_g = 2.0,
    double m_norm = 0.0,
    double switch_time = 0.0)
{
    // Список методов для сканирования
    struct MethodConfig {
        GuidanceMethod method;
        std::string name;
        double k;  // параметр k для метода
    };
    
    std::vector<MethodConfig> methods = {
        {GuidanceMethod::EVT, "EVT", 3.0},
        {GuidanceMethod::Gorka, "GORKA", 3.0}
        // Добавьте другие методы, если они есть
    };
    
    // Для каждого метода выполняем сканирование
    for (const auto& method : methods)
    {
        std::cout << "\n========================================" << std::endl;
        std::cout << "  СКАНИРОВАНИЕ ДЛЯ МЕТОДА: " << method.name << std::endl;
        std::cout << "========================================" << std::endl;
        
        // Диапазоны параметров
        const int MU_STEPS = 15;
        const int ETA_STEPS = 15;
        const int KM_STEPS = 8;
        const int KP_STEPS = 8;
        
        const double MU_MIN = 0.1;
        const double MU_MAX = 0.6;
        const double ETA_MIN = 10.0;
        const double ETA_MAX = 80.0;
        const double KM_MIN = 1.0;
        const double KM_MAX = 10.0;
        const double KP_MIN = 0.05;
        const double KP_MAX = 1.0;

        // Открываем файлы с учётом метода
        std::string prefix = method.name;
        std::ofstream heatmap_file(prefix + "_range_heatmap.dat");
        std::ofstream curves_file(prefix + "_range_curves.dat");
        std::ofstream full_data_file(prefix + "_range_full_data.dat");
        std::ofstream best_params_file(prefix + "_best_parameters.txt");
        
        if (!heatmap_file.is_open() || !curves_file.is_open() || 
            !full_data_file.is_open() || !best_params_file.is_open())
        {
            std::cerr << "❌ Не удалось создать файлы для записи" << std::endl;
            continue;
        }

        // Заголовки
        heatmap_file << "# Карта дальностей для метода " << method.name << std::endl;
        heatmap_file << "# Формат: mu eta max_range_km" << std::endl;
        
        curves_file << "# Кривые максимальной дальности для метода " << method.name << std::endl;
        curves_file << "# Формат: mu max_range_km eta_opt Km_opt Kp_opt" << std::endl;

        full_data_file << "# Полные данные сканирования для метода " << method.name << std::endl;
        full_data_file << "# Формат: mu eta Km Kp range_km" << std::endl;

        best_params_file << "ОПТИМАЛЬНЫЕ ПАРАМЕТРЫ ДЛЯ МЕТОДА " << method.name << std::endl;
        best_params_file << "================================================" << std::endl;

        std::cout << "\nПараметры сканирования для " << method.name << ":" << std::endl;
        std::cout << "  mu: [" << MU_MIN << ", " << MU_MAX << "] (" << MU_STEPS << " точек)" << std::endl;
        std::cout << "  eta: [" << ETA_MIN << ", " << ETA_MAX << "] (" << ETA_STEPS << " точек)" << std::endl;
        std::cout << "  K_m: [" << KM_MIN << ", " << KM_MAX << "] (" << KM_STEPS << " точек)" << std::endl;
        std::cout << "  K_P: [" << KP_MIN << ", " << KP_MAX << "] (" << KP_STEPS << " точек)" << std::endl;
        std::cout << "  K_g = " << K_g << " (фиксирован)" << std::endl;
        
        int total = MU_STEPS * ETA_STEPS * KM_STEPS * KP_STEPS;
        std::cout << "  Всего расчётов: " << total << std::endl;

        int completed = 0;
        int success_count = 0;

        // Массивы для хранения результатов
        std::map<std::pair<double, double>, double> max_range_for_mu_eta;
        std::map<std::pair<double, double>, std::tuple<double, double, double>> best_params_for_mu_eta;
        std::map<double, double> max_range_for_mu;
        std::map<double, std::tuple<double, double, double, double>> best_params_for_mu;

        // Для глобального максимума
        double global_max_range = 0.0;
        double global_mu = 0.0;
        double global_eta = 0.0;
        double global_Km = 0.0;
        double global_Kp = 0.0;
        TrajectoryResult global_best_result;

        auto start_time = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < MU_STEPS; ++i)
        {
            double mu = MU_MIN + i * (MU_MAX - MU_MIN) / (MU_STEPS - 1);
            
            for (int j = 0; j < ETA_STEPS; ++j)
            {
                double eta = ETA_MIN + j * (ETA_MAX - ETA_MIN) / (ETA_STEPS - 1);
                
                double best_range_for_pair = 0.0;
                double best_Km_for_pair = 0.0;
                double best_Kp_for_pair = 0.0;

                for (int k = 0; k < KM_STEPS; ++k)
                {
                    double Km = KM_MIN + k * (KM_MAX - KM_MIN) / (KM_STEPS - 1);
                    
                    for (int l = 0; l < KP_STEPS; ++l)
                    {
                        double Kp = KP_MIN + l * (KP_MAX - KP_MIN) / (KP_STEPS - 1);
                        
                        // Используем текущий метод
                        GuidanceMethod_Data methodData;
                        methodData.Method = method.method;
                        methodData.k = method.k;
                        methodData.K_g = K_g;

                        std::vector<TrajectoryParameters> s;
                        std::string resCode = CalcTrajectory_DZ(
                            calcParams, inits, limits,
                            CData, NData, RData,
                            methodData,
                            rocket,
                            mu, eta, Km, Kp, K_g,
                            s,
                            m_norm,
                            switch_time
                        );

                        completed++;
                        double range_km = 0.0;
                        
                        if (resCode == "0" && !s.empty())
                        {
                            double x_start = s.front().x_g_r;
                            double x_end = s.back().x_g_r;
                            range_km = (x_end - x_start) / 1000.0;
                            success_count++;
                        }

                        // Сохраняем полные данные
                        full_data_file << std::setprecision(6) 
                                       << mu << "\t" << eta << "\t" 
                                       << Km << "\t" << Kp << "\t" << range_km << "\n";

                        // Обновляем лучшую дальность для пары (mu, eta)
                        if (range_km > best_range_for_pair)
                        {
                            best_range_for_pair = range_km;
                            best_Km_for_pair = Km;
                            best_Kp_for_pair = Kp;
                        }

                        // Обновляем глобальный максимум
                        if (range_km > global_max_range)
                        {
                            global_max_range = range_km;
                            global_mu = mu;
                            global_eta = eta;
                            global_Km = Km;
                            global_Kp = Kp;
                            // Сохраняем результат
                            global_best_result.is_success = (resCode == "0" && !s.empty());
                            global_best_result.trajectory = s;
                            global_best_result.range_km = range_km;
                            global_best_result.flight_time = s.empty() ? 0.0 : s.back().t;
                            global_best_result.result_code = resCode;
                        }

                        // Прогресс
                        if (completed % 200 == 0 || completed == total)
                        {
                            auto current_time = std::chrono::high_resolution_clock::now();
                            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current_time - start_time).count();
                            
                            std::cout << "  Прогресс: " << completed << "/" << total 
                                      << " (" << std::fixed << std::setprecision(1) 
                                      << (double)completed / total * 100 << "%)"
                                      << ", успехов: " << success_count
                                      << ", время: " << elapsed << " с"
                                      << ", лучшая дальность: " << std::setprecision(2)
                                      << global_max_range << " км" << std::endl;
                        }
                    }
                }

                // Сохраняем лучшую дальность для пары (mu, eta)
                auto key = std::make_pair(mu, eta);
                max_range_for_mu_eta[key] = best_range_for_pair;
                best_params_for_mu_eta[key] = std::make_tuple(best_Km_for_pair, best_Kp_for_pair, best_range_for_pair);

                // Обновляем лучшую дальность для mu
                if (best_range_for_pair > max_range_for_mu[mu])
                {
                    max_range_for_mu[mu] = best_range_for_pair;
                    best_params_for_mu[mu] = std::make_tuple(eta, best_Km_for_pair, best_Kp_for_pair, best_range_for_pair);
                }
            }
        }

        auto end_time = std::chrono::high_resolution_clock::now();
        auto total_time = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();

        // Записываем карту дальностей
        for (const auto& item : max_range_for_mu_eta)
        {
            double mu = item.first.first;
            double eta = item.first.second;
            double range = item.second;
            heatmap_file << std::setprecision(6) << mu << "\t" << eta << "\t" << range << "\n";
        }

        // Записываем кривые максимальной дальности
        for (const auto& item : max_range_for_mu)
        {
            double mu = item.first;
            double range = item.second;
            auto params = best_params_for_mu[mu];
            double eta_opt = std::get<0>(params);
            double Km_opt = std::get<1>(params);
            double Kp_opt = std::get<2>(params);
            curves_file << std::setprecision(6) << mu << "\t" << range << "\t" 
                        << eta_opt << "\t" << Km_opt << "\t" << Kp_opt << "\n";
        }

        heatmap_file.close();
        curves_file.close();
        full_data_file.close();

        // ========================================================================
        // ВЫВОД РЕЗУЛЬТАТОВ ДЛЯ МЕТОДА
        // ========================================================================
        std::cout << "\n========================================" << std::endl;
        std::cout << "  РЕЗУЛЬТАТЫ ДЛЯ МЕТОДА: " << method.name << std::endl;
        std::cout << "========================================" << std::endl;
        std::cout << "Всего расчётов:      " << total << std::endl;
        std::cout << "Успешных встреч:     " << success_count << std::endl;
        std::cout << "Успешность:          " << std::fixed << std::setprecision(1) 
                  << (double)success_count / total * 100 << "%" << std::endl;
        std::cout << "Время выполнения:    " << total_time << " с" << std::endl;

        std::cout << "\n  ★ МАКСИМУМ ДАЛЬНОСТИ ДЛЯ " << method.name << " ★" << std::endl;
        std::cout << "  Дальность:          " << std::setprecision(2) << global_max_range << " км" << std::endl;
        std::cout << "  Время полёта:       " << global_best_result.flight_time << " с" << std::endl;
        std::cout << "\n  Оптимальные параметры:" << std::endl;
        std::cout << "    μ  = " << std::setprecision(4) << global_mu << std::endl;
        std::cout << "    η  = " << std::setprecision(2) << global_eta << std::endl;
        std::cout << "    K_m = " << std::setprecision(3) << global_Km << std::endl;
        std::cout << "    K_P = " << std::setprecision(4) << global_Kp << std::endl;
        std::cout << "    K_g = " << K_g << " (фиксирован)" << std::endl;

        // Запись в файл лучших параметров
        best_params_file << "\nМАКСИМУМ ДАЛЬНОСТИ" << std::endl;
        best_params_file << "==============================" << std::endl;
        best_params_file << "Дальность:          " << std::fixed << std::setprecision(2) 
                         << global_max_range << " км" << std::endl;
        best_params_file << "Время полёта:       " << global_best_result.flight_time << " с" << std::endl;
        best_params_file << "\nОптимальные параметры:" << std::endl;
        best_params_file << "  μ  = " << std::setprecision(4) << global_mu << std::endl;
        best_params_file << "  η  = " << std::setprecision(2) << global_eta << std::endl;
        best_params_file << "  K_m = " << std::setprecision(3) << global_Km << std::endl;
        best_params_file << "  K_P = " << std::setprecision(4) << global_Kp << std::endl;
        best_params_file << "  K_g = " << K_g << " (фиксирован)" << std::endl;
        best_params_file.close();

        // Сохраняем оптимальную траекторию для метода
        if (global_best_result.is_success && !global_best_result.trajectory.empty())
        {
            std::cout << "\n--- СОХРАНЕНИЕ ОПТИМАЛЬНОЙ ТРАЕКТОРИИ ДЛЯ " << method.name << " ---" << std::endl;
            
            std::string file_prefix = method.name + "_OPTIMAL_mu" + std::to_string(global_mu).substr(0, 5) +
                                      "_eta" + std::to_string(global_eta).substr(0, 4) +
                                      "_Km" + std::to_string(global_Km).substr(0, 4) +
                                      "_Kp" + std::to_string(global_Kp).substr(0, 5);
            
            saveAllResults(global_best_result.trajectory, file_prefix);
            saveAllPlots(global_best_result.trajectory, file_prefix,
                         method.name + " (дальность " +
                         std::to_string(global_max_range).substr(0, 5) + " км)",
                         inits.y_g_pusk);
            
            std::cout << "  ✅ Оптимальная траектория сохранена!" << std::endl;
        }
    }
}

// ============================================================================
// ГЛАВНАЯ ФУНКЦИЯ
// ============================================================================
int main(int argc, char* argv[])
{
    std::cout << "========================================" << std::endl;
    std::cout << "  РАСЧЁТ ТРАЕКТОРИИ РАКЕТЫ" << std::endl;
    std::cout << "  Вариант 26: Максимизация дальности" << std::endl;
    std::cout << "========================================" << std::endl;

    // ========================================================================
    // 1. Параметры расчёта
    // ========================================================================
    CalcParameters calcParams;
    calcParams.dt = 0.01;
    calcParams.r_por = 20.0;

    // ========================================================================
    // 2. Данные по цели, носителю и ракете
    // ========================================================================
    C_Data CData;
    CData.n_xa_c_potr = {{0.0, 0.0}, {600.0, 0.0}};
    CData.n_ya_c_potr = {{0.0, 0.0}, {600.0, 0.0}};

    N_Data NData;
    NData.n_xa_n_potr = {{0.0, 0.0}, {600.0, 0.0}};
    NData.n_ya_n_potr = {{0.0, 1.0}, {600.0, 1.0}};

    R_Data RData;
    RData.n_xa_r_potr = {{0.0, 0.0}, {600.0, 0.0}};

    // ========================================================================
    // 3. Параметры ракеты
    // ========================================================================
    RocketParams rocket;
    rocket.m_bch = 150.0;
    rocket.m_oun = 86.0;
    rocket.m_np = 48.0;
    rocket.d = 0.65;
    double m_t_init = 300.0;
    rocket.beta = calculate_beta(rocket.d, m_t_init);

    // ========================================================================
    // 4. Общие ограничения
    // ========================================================================
    LimitConditions limits;
    limits.n_ya_r_max = 40.0;
    limits.t_r_min = 2.0;
    limits.t_r_max = 600.0;
    limits.R_r_min = std::nan("");
    limits.R_r_max = std::nan("");
    limits.y_g_r_min = 0.0;
    limits.y_g_r_max = std::nan("");
    limits.epsilon_nc_min = std::nan("");
    limits.epsilon_nc_max = std::nan("");
    limits.v_r_min = 50.0;

    // ========================================================================
    // 5. НАЧАЛЬНЫЕ УСЛОВИЯ (ЦЕЛЬ НА ВСТРЕЧНОМ КУРСЕ)
    // ========================================================================
    InitConditions inits;
    inits.Theta_c0 = 180.0 * PI / 180.0;  // ВСТРЕЧНЫЙ КУРС
    inits.v_c0 = 18.0 * 0.514444;         // 18 узлов в м/с
    inits.x_g_c0 = 200000.0;              // 200 км от точки пуска
    inits.y_g_c0 = 0.0;
    inits.Theta_n0 = 0.0;
    inits.v_n0 = 120.0;
    inits.x_g_pusk = 0.0;
    inits.y_g_pusk = 200.0;
    inits.Theta_r0 = std::nan("");
    inits.v_r0 = 150.0;

    // ========================================================================
    // 6. СКАНИРОВАНИЕ И ОПТИМИЗАЦИЯ
    // ========================================================================
    scanAndOptimizeAllMethods(
        calcParams, inits, limits,
        CData, NData, RData,
        rocket,
        2.0,  // K_g фиксирован
        0.0,  // m_norm
        0.0   // switch_time
    );

    // ========================================================================
    // 7. ЗАВЕРШЕНИЕ
    // ========================================================================
    std::cout << "\n========================================" << std::endl;
    std::cout << "  РАСЧЁТ ЗАВЕРШЁН" << std::endl;
    std::cout << "========================================" << std::endl;

    std::cout << "\nДля построения графиков выполните:" << std::endl;
    std::cout << "  python plot_results.py" << std::endl;
    std::cout << "\nИли используйте Gnuplot:" << std::endl;
    std::cout << "  gnuplot -persist plots/OPTIMAL_*_plot.gnu" << std::endl;

    return 0;
}