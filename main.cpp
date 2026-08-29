#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <thread>
#include <mutex>
#include <atomic>

#include "CalcTrajectory.h"
#include "DataStructs.h"
#include "LibConstFunc.h"
#include "KinematicFunc2D.h"
#include "Atm_GOST4401.h"
#include "GuidanceMethods.h"

// ============================================================================
// СТРУКТУРА ДЛЯ ХРАНЕНИЯ УСПЕШНЫХ НАБОРОВ ПАРАМЕТРОВ
// ============================================================================
struct GorkaSuccessParams {
    double m0;
    double mu;
    double beta;
    double eta;
    double Km;
    double Kp;
    double Kg;
    double k;
    double delta_x_virt;
    double rocket_length;
    double t_gorka;
    bool passed_low_evt = false;
    bool passed_high_evt = false;
    double max_range = 0.0;
    double t_low = 0.0;
    double t_high = 0.0;
};

// ============================================================================
// ГЛОБАЛЬНЫЕ СИНХРОНИЗАЦИИ
// ============================================================================
std::mutex g_best_mutex;
std::mutex g_log_mutex;
std::mutex g_gorka_list_mutex;
std::vector<GorkaSuccessParams> g_gorka_success_list;

// ============================================================================
// РАБОЧИЙ ПОТОК: обрабатывает свой чанк mu
// ============================================================================
void worker_thread(
    size_t thread_id,
    size_t num_threads,
    const std::vector<double>& mu_vals,
    const std::vector<double>& eta_vals,
    const std::vector<double>& km_vals,
    const std::vector<double>& kp_vals,
    const std::vector<double>& kg_vals,
    const std::vector<double>& k_vals,
    const RocketParams& rocket_template,
    const CalcParameters& calc,
    const LimitConditions& limits,
    const C_Data& Cdata,
    const N_Data& Ndata,
    const R_Data& Rdata,
    const InitConditions& init_low,
    double range_min_km,
    double range_max_km,
    double range_step_km,
    OptimizeResult& global_best,
    std::atomic<int>& checked,
    std::atomic<int>& base_ok_count,
    std::atomic<int>& high_ok_count,
    std::ofstream& logFile)
{
    // Каждый поток берёт mu с шагом num_threads (чанкирование)
    for (size_t idx = thread_id; idx < mu_vals.size(); idx += num_threads)
    {
        double mu = mu_vals[idx];
        double m_pn = rocket_template.m_pn();
        double m_0 = m_pn / (1.0 - mu);
        double m_t = mu * m_0;
        double beta = 1.0;

        for (int iter = 0; iter < 5; ++iter)
        {
            beta = calculate_beta(rocket_template.d, m_t);
            if (beta * mu >= 1.0) break;
            m_0 = m_pn / (1.0 - beta * mu);
            m_t = mu * m_0;
        }

        if (beta * mu >= 1.0) continue;
        if (m_0 > 650.0) continue;

        const double ROCKET_DENSITY = 1555.0;      // кг/м³
        const double MAX_LENGTH = 4.1;             // м
        double rocket_volume = m_0 / ROCKET_DENSITY;
        double rocket_length = (24.0 * rocket_volume) / (5.0 * M_PI * rocket_template.d * rocket_template.d);

        if (rocket_length > MAX_LENGTH) {
            checked.fetch_add(1, std::memory_order_relaxed);
            continue;
        }

        RocketParams r = rocket_template;
        r.beta = beta;

        for (double eta : eta_vals)
        {
            for (double km : km_vals)
            {
                for (double kp : kp_vals)
                {
                    for (double kg : kg_vals)
                    {
                        for (double k_val : k_vals)
                        {
                            checked.fetch_add(1, std::memory_order_relaxed);

                            GuidanceMethod_Data mtd;
                            mtd.Method = GuidanceMethod::EVT;
                            mtd.k = k_val;
                            mtd.K_g = kg;

                            // === ШАГ 1: БАЗОВЫЙ СЦЕНАРИЙ ===
                            std::vector<TrajectoryParameters> s_low;
                            std::string code_low = CalcTrajectory(
                                calc, init_low, limits,
                                Cdata, Ndata, Rdata,
                                mtd, r,
                                mu, eta, km, kp, kg,
                                s_low);

                            if (code_low != "0" || s_low.empty() || s_low.back().t > limits.t_r_max)
                            {
                                continue;
                            }

                            base_ok_count.fetch_add(1, std::memory_order_relaxed);

                            // === ШАГ 2: ПЕРЕБОР ДАЛЬНОСТЕЙ ===
                            for (double range_km = range_min_km;
                                 range_km <= range_max_km;
                                 range_km += range_step_km)
                            {
                                double range_m = range_km * 1000.0;

                                InitConditions init_high;
                                init_high.x_g_pusk = 0.0;
                                init_high.y_g_pusk = 10000.0;
                                init_high.v_n0 = 270.0;
                                init_high.Theta_n0 = 0.0;
                                init_high.v_r0 = 270.0;
                                init_high.Theta_r0 = std::nan("");
                                init_high.x_g_c0 = range_m;
                                init_high.y_g_c0 = 0.0;
                                init_high.v_c0 = 18.0 * 0.514444;
                                init_high.Theta_c0 = M_PI;

                                std::vector<TrajectoryParameters> s_high;

                                std::string code_high = CalcTrajectory(
                                    calc, init_high, limits,
                                    Cdata, Ndata, Rdata,
                                    mtd, r,
                                    mu, eta, km, kp, kg,
                                    s_high);

                                if (code_high == "0" && !s_high.empty() && s_high.back().t <= limits.t_r_max)
                                {
                                    high_ok_count.fetch_add(1, std::memory_order_relaxed);

                                    // Потокобезопасное логирование
                                    {
                                        std::lock_guard<std::mutex> lock(g_log_mutex);
                                        logFile << "BOTH_OK: m0=" << m_0
                                                << " mu=" << mu << " beta=" << beta
                                                << " eta=" << eta << " Km=" << km
                                                << " Kp=" << kp << " Kg=" << kg << " k=" << k_val
                                                << " R_max=" << range_km
                                                << " t_low=" << s_low.back().t
                                                << " t_high=" << s_high.back().t << "\n";
                                        logFile.flush();
                                    }

                                    // Потокобезопасное обновление лучшего результата
                                    {
                                        std::lock_guard<std::mutex> lock(g_best_mutex);
                                        if (range_km > global_best.max_range)
                                        {
                                            global_best.max_range = range_km;
                                            global_best.m0 = m_0;
                                            global_best.mu = mu;
                                            global_best.beta = beta;
                                            global_best.eta = eta;
                                            global_best.Km = km;
                                            global_best.Kp = kp;
                                            global_best.Kg = kg;
                                            global_best.k = k_val;
                                            global_best.t_high = s_high.back().t;
                                            global_best.r_high = s_high.back().r_rc;
                                            global_best.t_low = s_low.back().t;
                                            global_best.r_low = s_low.back().r_rc;
                                            global_best.valid = true;
                                            global_best.rocket_length = rocket_length;

                                            std::cout << ">>> НОВЫЙ РЕКОРД: дальность=" << range_km
                                                      << " км, m0=" << m_0
                                                      << " mu=" << mu << " beta=" << beta
                                                      << " eta=" << eta << " Km=" << km
                                                      << " Kp=" << kp << " Kg=" << kg << " k=" << k_val
                                                      << " L=" << rocket_length << " м"
                                                      << " t_high=" << s_high.back().t
                                                      << " t_low=" << s_low.back().t << std::endl;
                                        }
                                    }
                                }
                            } // range_km

                            int current = checked.load(std::memory_order_relaxed);
                            if (current % 1000 == 0)
                            {
                                std::lock_guard<std::mutex> lock(g_log_mutex);
                                std::cout << "Прогресс: " << current
                                          << " | Базовых OK: " << base_ok_count.load()
                                          << " | Высоких OK: " << high_ok_count.load()
                                          << " | Рекорд: " << global_best.max_range << " км" << std::endl;
                            }

                        } // k_val
                    } // kg
                } // kp
            } // km
        } // eta
    } // mu (чанк)
}

// ============================================================================
// СТРУКТУРА ДЛЯ ХРАНЕНИЯ РЕЗУЛЬТАТОВ ОПТИМИЗАЦИИ ГОРКА
// ============================================================================
struct GorkaOptimizationResult {
    double range_km;
    double delta_x_virt;
    double t_flight;
    double v_impact;
    double r_final;
    bool success;
    std::string error_code;
};

// ============================================================================
// АЛГОРИТМ ОПТИМИЗАЦИИ ДЛЯ МАНЕВРА GORKA С ПЕРЕБОРОМ ДАЛЬНОСТЕЙ
// ============================================================================
std::vector<GorkaOptimizationResult> optimize_gorka_ranges(
    const RocketParams& rocket,
    const CalcParameters& calc,
    const LimitConditions& limits,
    const C_Data& Cdata,
    const N_Data& Ndata,
    const R_Data& Rdata,
    const InitConditions& init_Gorka,
    double mu, double eta, double km, double kp, double kg, double k,
    double range_start_km,
    double range_end_km,
    double range_step_km,
    const std::vector<double>& delta_x_virt_vals)
{
    std::vector<GorkaOptimizationResult> results;
    
    // Настройка метода Gorka с фиксированными параметрами
    GuidanceMethod_Data mtd_gorka;
    mtd_gorka.Method = GuidanceMethod::Gorka;
    mtd_gorka.K_H_march = 0.35;
    mtd_gorka.K_v_march = 0.38;
    mtd_gorka.K_H_gorka = 0.23;
    mtd_gorka.K_v_gorka = 0.40;
    mtd_gorka.H_march = 5.5;
    mtd_gorka.H_gorka = 200.0;
    mtd_gorka.r_save = 5000.0;
    mtd_gorka.k_dive = 4.0;
    mtd_gorka.H_save = 0.5;
    mtd_gorka.H_ant = 30.0;
    mtd_gorka.k = k;
    mtd_gorka.K_g = kg;
    
    std::cout << "\n========================================" << std::endl;
    std::cout << "=== ОПТИМИЗАЦИЯ МАНЕВРА GORKA ===" << std::endl;
    std::cout << "Фиксированные параметры:" << std::endl;
    std::cout << "  m0 = " << rocket.m0_start << " кг" << std::endl;
    std::cout << "  mu = " << mu << std::endl;
    std::cout << "  beta = " << rocket.beta << std::endl;
    std::cout << "  eta = " << eta << std::endl;
    std::cout << "  Km = " << km << std::endl;
    std::cout << "  Kp = " << kp << std::endl;
    std::cout << "  Kg = " << kg << std::endl;
    std::cout << "  k = " << k << std::endl;
    std::cout << "Диапазон дальностей: " << range_start_km << " - " << range_end_km << " км, шаг " << range_step_km << " км" << std::endl;
    std::cout << "Перебор delta_x_virt: ";
    for (double dx : delta_x_virt_vals) std::cout << dx/1000.0 << " ";
    std::cout << " км" << std::endl;
    std::cout << "========================================" << std::endl;
    
    // Открываем лог-файл
    std::ofstream logFile("gorka_optimization_v26.txt", std::ios::app);
    logFile << "\n=== GORKA OPTIMIZATION (FIXED PARAMS) ===\n";
    logFile << "m0=" << rocket.m0_start << ", mu=" << mu << ", beta=" << rocket.beta 
            << ", eta=" << eta << ", Km=" << km << ", Kp=" << kp 
            << ", Kg=" << kg << ", k=" << k << std::endl;
    logFile << "Range_km\tdelta_x_virt_km\tt_flight\tv_impact\tr_final\tStatus" << std::endl;
    logFile.flush();
    
    int total_combinations = static_cast<int>(
        ((range_end_km - range_start_km) / range_step_km + 1) * delta_x_virt_vals.size()
    );
    int checked = 0;
    int success_count = 0;
    
    std::cout << "\nНачало перебора... (всего " << total_combinations << " комбинаций)" << std::endl;
    
    for (double range_km = range_start_km; range_km <= range_end_km; range_km += range_step_km)
    {
        double range_m = range_km * 1000.0;
        
        // Обновляем начальные условия для текущей дальности
        InitConditions init_current = init_Gorka;
        init_current.x_g_c0 = range_m;
        
        for (double dx_val : delta_x_virt_vals)
        {
            checked++;
            
            // Устанавливаем delta_x_virt
            mtd_gorka.delta_x_virt = dx_val;
            
            std::vector<TrajectoryParameters> s_gorka;
            std::string code = CalcTrajectory(
                calc, init_current, limits,
                Cdata, Ndata, Rdata,
                mtd_gorka, rocket,
                mu, eta, km, kp, kg,
                s_gorka);
            
            GorkaOptimizationResult result;
            result.range_km = range_km;
            result.delta_x_virt = dx_val;
            result.success = false;
            result.error_code = code;
            
            if (code == "0" && !s_gorka.empty() && s_gorka.back().t <= limits.t_r_max)
            {
                result.success = true;
                result.t_flight = s_gorka.back().t;
                result.v_impact = s_gorka.back().v_r;
                result.r_final = s_gorka.back().r_rc;
                success_count++;
                
                // Вывод в консоль
                std::cout << "  R_Gorka = " << std::fixed << std::setprecision(1) << range_km 
                          << " км, delta_x_virt = " << std::setprecision(1) << dx_val/1000.0 
                          << " км, t = " << std::setprecision(2) << result.t_flight 
                          << " с, V = " << std::setprecision(1) << result.v_impact 
                          << " м/с, r_final = " << std::setprecision(1) << result.r_final << " м" << std::endl;
                
                // Запись в лог-файл
                logFile << std::fixed << std::setprecision(1) << range_km << "\t"
                        << std::setprecision(1) << dx_val/1000.0 << "\t"
                        << std::setprecision(2) << result.t_flight << "\t"
                        << std::setprecision(1) << result.v_impact << "\t"
                        << std::setprecision(1) << result.r_final << "\tOK\n";
                logFile.flush();
            }
            else
            {
                // Запись неудачных попыток в лог (опционально)
                logFile << std::fixed << std::setprecision(1) << range_km << "\t"
                        << std::setprecision(1) << dx_val/1000.0 << "\t"
                        << "-\t-\t-\tFAIL(" << code << ")\n";
            }
            
            results.push_back(result);
            
            // Прогресс
            if (checked % 10 == 0 || checked == total_combinations)
            {
                std::cout << "Прогресс: " << checked << "/" << total_combinations 
                          << " | Успешно: " << success_count << std::endl;
            }
        }
    }
    
    logFile << "\n=== DONE ===\n";
    logFile << "Total checked: " << checked << "\n";
    logFile << "Successful: " << success_count << "\n";
    logFile.close();
    
    // Вывод статистики
    std::cout << "\n========================================" << std::endl;
    std::cout << "=== СТАТИСТИКА ОПТИМИЗАЦИИ GORKA ===" << std::endl;
    std::cout << "Всего проверено комбинаций: " << checked << std::endl;
    std::cout << "Успешных: " << success_count << std::endl;
    
    // Находим максимальную дальность
    double max_range = 0.0;
    double max_dx = 0.0;
    double max_t = 0.0;
    double max_v = 0.0;
    
    for (const auto& r : results) {
        if (r.success && r.range_km > max_range) {
            max_range = r.range_km;
            max_dx = r.delta_x_virt;
            max_t = r.t_flight;
            max_v = r.v_impact;
        }
    }
    
    if (max_range > 0) {
        std::cout << "\n=== ЛУЧШИЙ РЕЗУЛЬТАТ ===" << std::endl;
        std::cout << "Максимальная дальность Gorka: " << max_range << " км" << std::endl;
        std::cout << "  delta_x_virt = " << max_dx/1000.0 << " км" << std::endl;
        std::cout << "  Время полета = " << max_t << " с" << std::endl;
        std::cout << "  Скорость подлета = " << max_v << " м/с" << std::endl;
    }
    std::cout << "========================================" << std::endl;
    
    return results;
}


// ============================================================================
// ФУНКЦИЯ ВЫВОДА РЕЗУЛЬТАТОВ ОПТИМИЗАЦИИ GORKA В ТАБЛИЦЕ
// ============================================================================
void print_gorka_results_table(const std::vector<GorkaOptimizationResult>& results) {
    if (results.empty()) {
        std::cout << "\nНет результатов для отображения." << std::endl;
        return;
    }
    
    std::cout << "\n=== ТАБЛИЦА РЕЗУЛЬТАТОВ GORKA ===" << std::endl;
    std::cout << std::setw(12) << "R, км" 
              << std::setw(18) << "dx_virt, км"
              << std::setw(12) << "t, с"
              << std::setw(12) << "V, м/с"
              << std::setw(12) << "r_final, м"
              << std::setw(10) << "Status" << std::endl;
    std::cout << std::string(76, '-') << std::endl;
    
    // Сортируем по дальности
    auto sorted = results;
    std::sort(sorted.begin(), sorted.end(),
        [](const GorkaOptimizationResult& a, const GorkaOptimizationResult& b) {
            return a.range_km < b.range_km;
        });
    
    for (const auto& r : sorted) {
        std::cout << std::setw(12) << std::fixed << std::setprecision(1) << r.range_km
                  << std::setw(18) << std::setprecision(1) << r.delta_x_virt/1000.0;
        
        if (r.success) {
            std::cout << std::setw(12) << std::setprecision(2) << r.t_flight
                      << std::setw(12) << std::setprecision(1) << r.v_impact
                      << std::setw(12) << std::setprecision(1) << r.r_final
                      << std::setw(10) << "OK";
        } else {
            std::cout << std::setw(12) << "-"
                      << std::setw(12) << "-"
                      << std::setw(12) << "-"
                      << std::setw(10) << r.error_code;
        }
        std::cout << std::endl;
    }
}

// ============================================================================
// ОСНОВНАЯ ФУНКЦИЯ ОПТИМИЗАЦИИ (распараллеленная)
// ============================================================================
OptimizeResult optimize_max_range_parallel(
    const RocketParams& rocket_template,
    const CalcParameters& calc,
    const LimitConditions& limits,
    const C_Data& Cdata,
    const N_Data& Ndata,
    const R_Data& Rdata,
    const InitConditions& init_low,
    double range_min_km,
    double range_max_km,
    double range_step_km)
{
    OptimizeResult global_best;
    global_best.valid = false;
    global_best.max_range = 0.0;

    // === СЕТКА ПАРАМЕТРОВ ===
    std::vector<double> mu_vals;
    for (double mu = 0.44; mu <= 0.6501; mu += 0.01)
        mu_vals.push_back(mu);

    std::vector<double> eta_vals = {6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0, 24.0};
    std::vector<double> km_vals = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> kp_vals = {0.25, 0.30, 0.35, 0.40, 0.45, 0.50};
    std::vector<double> kg_vals = {1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
    std::vector<double> k_vals = {1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0};

    int total_combinations = static_cast<int>(
        mu_vals.size() * eta_vals.size() * km_vals.size() *
        kp_vals.size() * kg_vals.size() * k_vals.size()
    );

    std::atomic<int> checked{0};
    std::atomic<int> base_ok_count{0};
    std::atomic<int> high_ok_count{0};

    std::ofstream logFile("optimization_v26.txt", std::ios::app);
    logFile << "\n=== V26 OPTIMIZATION (MAX RANGE) PARALLEL START ===\n";
    logFile << "Threads: 16\n";
    logFile << "Total param combinations: " << total_combinations << "\n";
    logFile << "Base scenario: 0.2 km, 120 m/s, 150 km\n";
    logFile << "Range scan: " << range_min_km << "-" << range_max_km << " km, step " << range_step_km << "\n";
    logFile << "Mass limit: m0 <= 650 kg\n";
    logFile << "Rocket length limit: <= 4.1 m (density 1555 kg/m³)\n";
    logFile << "Grid: mu[0.44:0.02:0.70] (outer loop), eta[6:18], Km[1:5], Kp[0.25:0.05:0.80], Kg[1.05:0.1:2.0], k[3:10]\n";
    logFile.flush();

    std::cout << "=== НАЧАЛО ОПТИМИЗАЦИИ (16 потоков) ===" << std::endl;
    std::cout << "Всего комбинаций параметров: " << total_combinations << std::endl;
    std::cout << "Базовый сценарий: 0.2 km, 120 m/s, 150 km" << std::endl;
    std::cout << "Перебор дальностей: " << range_min_km << "-" << range_max_km << " км" << std::endl;
    std::cout << "Ограничение по массе: m0 <= 650 кг" << std::endl;
    std::cout << "μ — внешний цикл (чанкирование по потокам)" << std::endl;
    std::cout << std::endl;

    const size_t num_threads = 16;
    std::vector<std::thread> threads;

    for (size_t t = 0; t < num_threads; ++t)
    {
        threads.emplace_back(worker_thread,
                             t, num_threads,
                             std::cref(mu_vals),
                             std::cref(eta_vals),
                             std::cref(km_vals),
                             std::cref(kp_vals),
                             std::cref(kg_vals),
                             std::cref(k_vals),
                             std::cref(rocket_template),
                             std::cref(calc),
                             std::cref(limits),
                             std::cref(Cdata),
                             std::cref(Ndata),
                             std::cref(Rdata),
                             std::cref(init_low),
                             range_min_km, range_max_km, range_step_km,
                             std::ref(global_best),
                             std::ref(checked),
                             std::ref(base_ok_count),
                             std::ref(high_ok_count),
                             std::ref(logFile));
    }

    for (auto& t : threads)
        t.join();

    logFile << "\n=== DONE ===\n";
    logFile << "Checked: " << checked.load() << "/" << total_combinations << "\n";
    logFile << "Base OK: " << base_ok_count.load() << " High OK: " << high_ok_count.load() << "\n";
    logFile.close();

    std::cout << "\n========================================" << std::endl;
    std::cout << "ОПТИМИЗАЦИЯ ЗАВЕРШЕНА" << std::endl;
    std::cout << "Проверено комбинаций: " << checked.load() << "/" << total_combinations << std::endl;
    std::cout << "Успешных базовых: " << base_ok_count.load() << std::endl;
    std::cout << "Успешных высоких: " << high_ok_count.load() << std::endl;

    if (global_best.valid)
    {
        std::cout << "\n=== ЛУЧШИЙ РЕЗУЛЬТАТ ===" << std::endl;
        std::cout << "Максимальная дальность: " << global_best.max_range << " км" << std::endl;
        std::cout << "m0 = " << global_best.m0 << " кг" << std::endl;
        std::cout << "mu = " << global_best.mu << std::endl;
        std::cout << "beta = " << global_best.beta << std::endl;
        std::cout << "eta = " << global_best.eta << std::endl;
        std::cout << "Km = " << global_best.Km << std::endl;
        std::cout << "Kp = " << global_best.Kp << std::endl;
        std::cout << "Kg = " << global_best.Kg << std::endl;
        std::cout << "k = " << global_best.k << std::endl;
        std::cout << "t_low = " << global_best.t_low << " с" << std::endl;
        std::cout << "t_high = " << global_best.t_high << " с" << std::endl;
        std::cout << "Длина ракеты = " << global_best.rocket_length << " м" << std::endl;
    }
    else
    {
        std::cout << "\nРешение не найдено!" << std::endl;
    }
    std::cout << "========================================" << std::endl;

    return global_best;
}
// ============================================================================
// ФУНКЦИЯ ПОСТ-ОБРАБОТКИ С СОЗДАНИЕМ ФАЙЛА РАДИОГОРИЗОНТА
// ============================================================================
void run_and_save(
    const std::string& case_name,
    const CalcParameters& calc,
    const LimitConditions& limits,
    const C_Data& Cdata,
    const N_Data& Ndata,
    const R_Data& Rdata,
    const InitConditions& init,
    GuidanceMethod_Data method,
    const RocketParams& rocket,
    double mu, double eta, double km, double kp, double kg)
{
    std::vector<TrajectoryParameters> result;
    std::string code = CalcTrajectory(
        calc, init, limits,
        Cdata, Ndata, Rdata,
        method, rocket,
        mu, eta, km, kp, kg,
        result);

    std::cout << "[" << case_name << "] Результат: " << code << std::endl;
    if (!result.empty())
    {
        std::cout << "  Время полета: " << result.back().t << " с" << std::endl;
        std::cout << "  Конечная дальность: " << result.back().x_g_r << " м" << std::endl;
        std::cout << "  Конечная высота: " << result.back().y_g_r << " м" << std::endl;
        std::cout << "  Конечная скорость: " << result.back().v_r << " м/с" << std::endl;
    }

    // === Сохранение в файлы ===
    std::ofstream C_File(case_name + "_C.gra");
    std::ofstream N_File(case_name + "_N.gra");
    std::ofstream R_File(case_name + "_R.gra");
    std::ofstream Kinematic_File(case_name + "_Kinematic.gra");
    std::ofstream Horizon_File(case_name + "_Horizon.gra"); // Файл радиогоризонта

    // Заголовки
    C_File << "t\tTheta_c\tv_c\tx_g_c\ty_g_c\tn_xa_c\tn_ya_c";
    N_File << "t\tTheta_n\tv_n\tx_g_n\ty_g_n\tn_xa_n\tn_ya_n";
    R_File << "t\tTheta_r\tv_r\tx_g_r\ty_g_r\tn_xa_r\tn_ya_r\tmass_r\tP_r\talpha_r";
    Kinematic_File << "t\tr_rc\tdot_r_rc\tepsilon_rc\tdot_epsilon_rc"
                   << "\tr_nr\tdot_r_nr\tepsilon_nr\tdot_epsilon_nr"
                   << "\tr_nc\tdot_r_nc\tepsilon_nc\tdot_epsilon_nc";
    Horizon_File << "t\ty_g_r\tr_rc\ty_radar_horizon\tis_visible";

    for (size_t i = 0; i < result.size(); i++)
    {
        N_File << "\n" << result[i].t << "\t" << result[i].Theta_n << "\t"
               << result[i].v_n << "\t" << result[i].x_g_n << "\t"
               << result[i].y_g_n << "\t" << result[i].n_xa_n << "\t"
               << result[i].n_ya_n;

        R_File << "\n" << result[i].t << "\t" << result[i].Theta_r << "\t"
               << result[i].v_r << "\t" << result[i].x_g_r << "\t"
               << result[i].y_g_r << "\t" << result[i].n_xa_r << "\t"
               << result[i].n_ya_r << "\t" << result[i].mass_r << "\t"
               << result[i].P_r << "\t" << result[i].alpha_r;

        C_File << "\n" << result[i].t << "\t" << result[i].Theta_c << "\t"
               << result[i].v_c << "\t" << result[i].x_g_c << "\t"
               << result[i].y_g_c << "\t" << result[i].n_xa_c << "\t"
               << result[i].n_ya_c;

        Kinematic_File << "\n" << result[i].t << "\t"
                       << result[i].r_rc << "\t" << result[i].dot_r_rc << "\t"
                       << result[i].epsilon_rc << "\t" << result[i].dot_epsilon_rc << "\t"
                       << result[i].r_nr << "\t" << result[i].dot_r_nr << "\t"
                       << result[i].epsilon_nr << "\t" << result[i].dot_epsilon_nr << "\t"
                       << result[i].r_nc << "\t" << result[i].dot_r_nc << "\t"
                       << result[i].epsilon_nc << "\t" << result[i].dot_epsilon_nc;

        // Запись данных радиогоризонта
        Horizon_File << "\n" << result[i].t << "\t"
                     << result[i].y_g_r << "\t"
                     << result[i].r_rc << "\t";
        
        if (result[i].y_radar_horizon.has_value()) {
            Horizon_File << result[i].y_radar_horizon.value() << "\t1";
        } else {
            Horizon_File << "-1\t0";
        }
    }

    C_File.close();
    N_File.close();
    R_File.close();
    Kinematic_File.close();
    Horizon_File.close();

    std::cout << "  Файлы сохранены: " << case_name << "_*.gra" << std::endl;
    std::cout << "  Файл радиогоризонта: " << case_name << "_Horizon.gra" << std::endl;
}

// ============================================================================
// MAIN
// ============================================================================
int main()
{
    C_Data Cdata;
    N_Data Ndata;
    R_Data Rdata;

    Cdata.n_xa_c_potr = {{0.0, 0.0}, {10.0, 0.0}};
    Cdata.n_ya_c_potr = {{0.0, -1.0}, {10.0, -1.0}};

    Ndata.n_xa_n_potr = {{0.0, 0.0}, {10.0, 0.0}};
    Ndata.n_ya_n_potr = {{0.0, 1.0}, {10.0, 1.0}};

    Rdata.n_xa_r_potr = {{0.0, 0.0}, {10.0, 0.0}};

    RocketParams rocket;
    rocket.m_bch = 150.0;
    rocket.m_oun = 65.0;
    rocket.m_np  = 11.0;
    rocket.d     = 0.400;
    rocket.beta  = 1.13773; // ИЗМЕНИТЬ
    
    // Расчет m0_start для корректной работы
    double m_pn = rocket.m_pn();
    double mu_fixed = 0.57;
    rocket.m0_start = m_pn / (1.0 - rocket.beta * mu_fixed);

    CalcParameters calc;
    calc.dt = 0.01;
    calc.r_por = 15.0;

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

    InitConditions init_low;
    init_low.Theta_c0 = M_PI;
    init_low.v_c0 = 18.0 * 0.514444;
    init_low.x_g_c0 = 300000.0; // ИЗМЕНИТЬ
    init_low.y_g_c0 = 0.0;
    init_low.Theta_n0 = 0.0;
    init_low.v_n0 = 120.0;
    init_low.x_g_pusk = 0.0;
    init_low.y_g_pusk = 200.0;
    init_low.Theta_r0 = std::nan("");
    init_low.v_r0 = init_low.v_n0;

    InitConditions init_high;
    init_high.x_g_pusk = 0.0;
    init_high.y_g_pusk = 10000.0; 
    init_high.v_n0 = 270.0;
    init_high.Theta_n0 = 0.0;
    init_high.v_r0 = 270.0;
    init_high.Theta_r0 = std::nan("");
    init_high.x_g_c0 = 1040000.0; // ИЗМЕНИТЬ
    init_high.y_g_c0 = 0.0;
    init_high.v_c0 = 18.0 * 0.514444;
    init_high.Theta_c0 = M_PI;

    InitConditions init_Gorka;
    init_Gorka.Theta_c0 = M_PI;
    init_Gorka.v_c0 = 18.0 * 0.514444;
    init_Gorka.x_g_c0 =50000.0; // ИЗМЕНИТЬ
    init_Gorka.y_g_c0 = 0.0;
    init_Gorka.Theta_n0 = 0.0;
    init_Gorka.v_n0 = 120.0;
    init_Gorka.x_g_pusk = 0.0;
    init_Gorka.y_g_pusk = 200.0;
    init_Gorka.Theta_r0 = std::nan("");
    init_Gorka.v_r0 = init_Gorka.v_n0;

    GuidanceMethod_Data method;
    method.Method = GuidanceMethod::EVT;
    method.k = 5.0; // ИЗМЕНИТЬ
    method.K_g = 1.5; // ИЗМЕНИТЬ

    // === 1. ЗАПУСК ОПТИМИЗАЦИИ (закомментирован, так как параметры уже найдены) ===
    
    // OptimizeResult result = optimize_max_range_parallel(
    //     rocket, calc, limits,
    //     Cdata, Ndata, Rdata,
    //     init_low,
    //     1000.0,
    //     1150.0,
    //     10.0
    // );


    // BOTH_OK: m0=515.099 mu=0.49 beta=1.14541 eta=6 Km=3 Kp=0.25 Kg=1.9 k=3.5 R_max=700 t_low=277.35 t_high=544.93
    // BOTH_OK: m0=515.099 mu=0.49 beta=1.14541 eta=6 Km=5 Kp=0.25 Kg=2 k=5 R_max=800 t_low=394.68 t_high=507.4 +++
    // m0=642.965 mu=0.57 beta=1.13773 eta=6 Km=5 Kp=0.25 Kg=1.5 k=5 R_max=1040 t_low=460.51 t_high=596.33

    // === ПАРАМЕТРЫ ===
    double mu = 0.57; // ИЗМЕНИТЬ
    double eta = 6.0; // ИЗМЕНИТЬ
    double Km = 5.0; // ИЗМЕНИТЬ
    double Kp = 0.25; // ИЗМЕНИТЬ
    double Kg = 1.5; // ИЗМЕНИТЬ
    double k = 5.0; // ИЗМЕНИТЬ

    // ===  ОПТИМИЗАЦИЯ GORKA С ПЕРЕБОРОМ ДАЛЬНОСТЕЙ ===

    std::vector<double> delta_x_virt_vals = {
        10000.0, 15000.0, 20000.0, 25000.0, 30000.0, 35000.0, 40000.0, 45000.0, 
        50000.0, 55000.0, 60000.0, 65000.0, 70000.0, 75000.0, 80000.0, 85000.0, 90000.0
    };

    // std::vector<GorkaOptimizationResult> gorka_results = optimize_gorka_ranges(
    //     rocket, calc, limits,
    //     Cdata, Ndata, Rdata,
    //     init_Gorka,
    //     mu, eta, Km, Kp, Kg, k,
    //     20.0,
    //     100.0,
    //     10,
    //     delta_x_virt_vals
    // );

    // // === ВЫВОД ТАБЛИЦЫ РЕЗУЛЬТАТОВ ===
    // print_gorka_results_table(gorka_results);

    // // === СОХРАНЕНИЕ ТРАЕКТОРИЙ ДЛЯ ЛУЧШИХ РЕЗУЛЬТАТОВ ===
    // // Находим лучший результат
    // double best_range = 0.0;
    // double best_dx = 0.0;
    // for (const auto& r : gorka_results) {
    //     if (r.success && r.range_km > best_range) {
    //         best_range = r.range_km;
    //         best_dx = r.delta_x_virt;
    //     }
    // }

    // === 1. НИЗКОВЫСОТНЫЙ СЦЕНАРИЙ (EVT) ===
    method.Method = GuidanceMethod::EVT;
    run_and_save("V26_Best_low", calc, limits,
                 Cdata, Ndata, Rdata, init_low,
                 method, rocket,
                 mu, eta, Km, Kp, Kg);

    // === 2. ВЫСОКОВЫСОТНЫЙ СЦЕНАРИЙ (EVT) ===
    method.Method = GuidanceMethod::EVT;
    run_and_save("V26_Best_high", calc, limits,
                 Cdata, Ndata, Rdata, init_high,
                 method, rocket,
                 mu, eta, Km, Kp, Kg);

    // === 3. МЕТОД ГОРКИ (GORKA) ===
    method.Method = GuidanceMethod::Gorka;
    method.K_H_march = 0.35;
    method.K_v_march = 0.38;
    method.K_H_gorka = 0.23;
    method.K_v_gorka = 0.40;
    method.H_march = 5.5;
    method.H_gorka = 200.0;
    method.r_save = 5000.0;
    method.k_dive = 4.0;
    method.H_save = 0.5;
    method.H_ant = 30.0;
    method.delta_x_virt = 10000.0; // ИЗМЕНИТЬ
    method.k = 5.0;
    method.K_g = 1.9;

    run_and_save("V26_Best_gorka", calc, limits,
                 Cdata, Ndata, Rdata, init_Gorka,
                 method, rocket,
                 mu, eta, Km, Kp, Kg);

    return 0;
}