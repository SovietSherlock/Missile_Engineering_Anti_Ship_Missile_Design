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
// ГЛОБАЛЬНЫЕ СИНХРОНИЗАЦИИ
// ============================================================================
std::mutex g_best_mutex;
std::mutex g_log_mutex;

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
        double rocket_length = (8.0 * rocket_volume) / (M_PI * rocket_template.d * rocket_template.d);
        
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
                            double mn_l = 0, st_l = 0;

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
    // μ: 0.44 до 0.70 с шагом 0.02 — САМЫЙ ВЕРХНИЙ ЦИКЛ
    std::vector<double> mu_vals;
    for (double mu = 0.44; mu <= 0.4901; mu += 0.01)
        mu_vals.push_back(mu);

    std::vector<double> eta_vals = {6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
    std::vector<double> km_vals = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> kp_vals = {0.25, 0.30, 0.35, 0.40, 0.45, 0.50};
    std::vector<double> kg_vals = {1.05, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
    std::vector<double> k_vals = {3, 4, 5, 6, 7, 8, 9, 10};

    int total_combinations = static_cast<int>(
        mu_vals.size() * eta_vals.size() * km_vals.size() *
        kp_vals.size() * kg_vals.size() * k_vals.size()
        );

    std::atomic<int> checked{0};
    std::atomic<int> base_ok_count{0};
    std::atomic<int> high_ok_count{0};

    // === ЛОГ-ФАЙЛ ===
    std::ofstream logFile("optimization_v26.txt", std::ios::app);
    logFile << "\n=== V26 OPTIMIZATION (MAX RANGE) PARALLEL START ===\n";
    logFile << "Threads: 16\n";
    logFile << "Total param combinations: " << total_combinations << "\n";
    logFile << "Base scenario: 0.2 km, 120 m/s, 150 km\n";
    logFile << "Range scan: " << range_min_km << "-" << range_max_km << " km, step " << range_step_km << "\n";
    logFile << "Mass limit: m0 <= 650 kg\n";
    logFile << "Rocket length limit: <= 4.1 m (density 1555 kg/m³)\n"; // <-- ДОБАВЛЕНО
    logFile << "Grid: mu[0.44:0.02:0.70] (outer loop), eta[6:18], Km[1:5], Kp[0.25:0.05:0.80], Kg[1.05:0.1:2.0], k[3:10]\n";
    logFile.flush();

    std::cout << "=== НАЧАЛО ОПТИМИЗАЦИИ (16 потоков) ===" << std::endl;
    std::cout << "Всего комбинаций параметров: " << total_combinations << std::endl;
    std::cout << "Базовый сценарий: 0.2 km, 120 m/s, 150 km" << std::endl;
    std::cout << "Перебор дальностей: " << range_min_km << "-" << range_max_km << " км" << std::endl;
    std::cout << "Ограничение по массе: m0 <= 650 кг" << std::endl;
    std::cout << "μ — внешний цикл (чанкирование по потокам)" << std::endl;
    std::cout << std::endl;

    // === ЗАПУСК 16 ПОТОКОВ ===
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

    // === ОЖИДАНИЕ ЗАВЕРШЕНИЯ ===
    for (auto& t : threads)
        t.join();

    // === ИТОГОВЫЙ ВЫВОД В ЛОГ ===
    logFile << "\n=== DONE ===\n";
    logFile << "Checked: " << checked.load() << "/" << total_combinations << "\n";
    logFile << "Base OK: " << base_ok_count.load() << " High OK: " << high_ok_count.load() << "\n";

    if (global_best.valid)
    {
        logFile << "\nBEST RESULT:\n";
        logFile << "Max range = " << global_best.max_range << " km\n";
        logFile << "m0 = " << global_best.m0 << " kg\n";
        logFile << "mu = " << global_best.mu << "\n";
        logFile << "beta = " << global_best.beta << "\n";
        logFile << "eta = " << global_best.eta << "\n";
        logFile << "Km = " << global_best.Km << "\n";
        logFile << "Kp = " << global_best.Kp << "\n";
        logFile << "Kg = " << global_best.Kg << "\n";
        logFile << "k = " << global_best.k << "\n";
        logFile << "t_low = " << global_best.t_low << " s\n";
        logFile << "t_high = " << global_best.t_high << " s\n";
        logFile << "Rocket length = " << global_best.rocket_length << " m\n"; 
    }
    else
    {
        logFile << "\nNo solution found!\n";
    }
    logFile << "========================================\n";
    logFile.flush();
    logFile.close();

    // === ИТОГОВЫЙ ВЫВОД В КОНСОЛЬ ===
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
        std::cout << "Длина ракеты = " << global_best.rocket_length << " м" << std::endl; // <-- ДОБАВЛЕНО
    }
    else
    {
        std::cout << "\nРешение не найдено!" << std::endl;
    }
    std::cout << "========================================" << std::endl;

    return global_best;
}


void run_and_save(
    const std::string& case_name,
    const CalcParameters& calc,
    const LimitConditions& limits,
    const C_Data& Cdata, const N_Data& Ndata, const R_Data& Rdata,
    const InitConditions& init,
    GuidanceMethod_Data method,
    const RocketParams& rocket,
    double mu, double eta, double km, double kp, double kg)
{
    std::vector<TrajectoryParameters> result;
    std::string code = CalcTrajectory(calc, init, limits, Cdata, Ndata, Rdata,
                                      method, rocket, mu, eta, km, kp, kg, result);

    std::cout << "[" << case_name << "] Результат: " << code << std::endl;
    if (!result.empty())
    {
        std::cout << "  Время полета: " << result.back().t << " с" << std::endl;
        std::cout << "  Конечная дальность: " << result.back().x_g_r << " м" << std::endl;
        std::cout << "  Конечная высота: " << result.back().y_g_r << " м" << std::endl;
        std::cout << "  Конечная скорость: " << result.back().v_r << " м/с" << std::endl;
    }

    std::ofstream C_File(case_name + "_C.gra");
    std::ofstream N_File(case_name + "_N.gra");
    std::ofstream R_File(case_name + "_R.gra");
    std::ofstream Kinematic_File(case_name + "_Kinematic.gra");

    C_File << "t\tTheta_c\tv_c\tx_g_c\ty_g_c\tn_xa_c\tn_ya_c";
    N_File << "t\tTheta_n\tv_n\tx_g_n\ty_g_n\tn_xa_n\tn_ya_n";
    R_File << "t\tTheta_r\tv_r\tx_g_r\ty_g_r\tn_xa_r\tn_ya_r\tmass_r\tP_r\talpha_r";
    Kinematic_File << "t\tr_rc\tdot_r_rc\tepsilon_rc\tdot_epsilon_rc"
                   << "\tr_nr\tdot_r_nr\tepsilon_nr\tdot_epsilon_nr"
                   << "\tr_nc\tdot_r_nc\tepsilon_nc\tdot_epsilon_nc";

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
    }

    N_File.close(); C_File.close(); R_File.close(); Kinematic_File.close();
    std::cout << "  Файлы сохранены: " << case_name << "_*.gra" << std::endl;
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
    rocket.m_oun = 86.0;
    rocket.m_np  = 48.0;
    rocket.d     = 0.510;
    rocket.beta  = 1.15763; //ИЗМЕНИШЬ ЗДЕСЬ БЕТУ НА НАЙДЕННУЮ

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
    init_low.Theta_c0 =  M_PI ;
    init_low.v_c0 = 18.0 * 0.514444;
    init_low.x_g_c0 = 150000.0;
    init_low.y_g_c0 = 0.0;
    init_low.Theta_n0 = 0.0;
    init_low.v_n0 = 120.0;              // НЕБЛАГОПРИЯТНЫЕ УСЛОВИЯ
    init_low.x_g_pusk = 0.0;
    init_low.y_g_pusk = 200.0;        // ВЫСОТА 200 м
    init_low.Theta_r0 = std::nan("");
    init_low.v_r0 = init_low.v_n0;


    InitConditions init_high;
    init_high.x_g_pusk = 0.0;
    init_high.y_g_pusk = 10000.0;
    init_high.v_n0 = 270.0;
    init_high.Theta_n0 = 0.0;
    init_high.v_r0 = 270.0;
    init_high.Theta_r0 = std::nan("");
    init_high.x_g_c0 = 1210000; // ЗДЕСБ ИЗМЕНИШЬ ДАЛЬНОСТЬ НА НАЙДЕННУЮ
    init_high.y_g_c0 = 0.0;
    init_high.v_c0 = 18.0 * 0.514444;
    init_high.Theta_c0 = M_PI;


    // OptimizeResult result = optimize_max_range_parallel(
    //     rocket, calc, limits,
    //     Cdata, Ndata, Rdata,
    //     init_low,
    //     1150.0,
    //     1300.0,
    //     10.0
    //     );

    double mu = 0.48, beta = 1.15763, eta = 6, Km = 5, Kp = 0.25, Kg = 1.9; // изменить при получении оптимизированных данных
    GuidanceMethod_Data method;
    method.Method = GuidanceMethod::EVT;
    method.k = 9;

    std::cout << "=== init_low ===" << std::endl;
    std::cout << "x_g_pusk=" << init_low.x_g_pusk << " y_g_pusk=" << init_low.y_g_pusk << std::endl;
    std::cout << "v_n0=" << init_low.v_n0 << " Theta_n0=" << init_low.Theta_n0 << std::endl;
    std::cout << "v_r0=" << init_low.v_r0 << " Theta_r0=" << init_low.Theta_r0 << std::endl;
    std::cout << "x_g_c0=" << init_low.x_g_c0 << " y_g_c0=" << init_low.y_g_c0 << std::endl;
    std::cout << "v_c0=" << init_low.v_c0 << " Theta_c0=" << init_low.Theta_c0 << std::endl;
    std::cout << "=== limits ===" << std::endl;
    std::cout << "t_r_max=" << limits.t_r_max << " n_ya_r_max=" << limits.n_ya_r_max << std::endl;
    std::cout << "=== calc ===" << std::endl;
    std::cout << "dt=" << calc.dt << " r_por=" << calc.r_por << std::endl;
    std::cout << "=== rocket ===" << std::endl;
    std::cout << "m_bch=" << rocket.m_bch << " m_oun=" << rocket.m_oun << " m_np=" << rocket.m_np << std::endl;
    std::cout << "d=" << rocket.d << " beta=" << rocket.beta << std::endl;
    std::cout << "m_pn=" << rocket.m_pn() << std::endl;
    std::cout << "=== method ===" << std::endl;
    std::cout << "k=" << method.k << std::endl;

    // Вычислите m_0 так же, как в оптимизаторе
    double m_pn = rocket.m_pn();
    double m_0_check = m_pn / (1.0 - beta * mu);
    double m_t = mu * m_0_check;
    for (int iter = 0; iter < 5; ++iter)
    {
        beta = calculate_beta(rocket.d, m_t);
        if (beta * mu >= 1.0) break;
        m_0_check = m_pn / (1.0 - beta * mu);
        m_t = mu * m_0_check;
    }
    std::cout << std::setprecision(17) << "Final beta = " << beta << std::endl;
    std::cout << std::setprecision(17) << "Final mu = " << mu << std::endl;
    std::cout << std::setprecision(17) << "Final mu*beta = " << mu * beta << std::endl;
    std::cout << std::setprecision(17) << "m_0_check = " << m_0_check << std::endl;
    std::cout << std::setprecision(17) << "m_0 from log = 639.15" << std::endl;
    rocket.beta = beta;
    rocket.beta = beta;

    run_and_save("V12_Best_low", calc, limits,
                 Cdata, Ndata, Rdata, init_low,
                 method, rocket, mu, eta, Km, Kp, Kg);

    run_and_save("V12_Best_high", calc, limits,
                 Cdata, Ndata, Rdata, init_high,
                 method, rocket, mu, eta, Km, Kp, Kg);

    return 0;




}
