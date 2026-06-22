// CalcTrajectory.cpp
#include "CalcTrajectory.h"
#include "KinematicFunc2D.h"
#include "GuidanceMethods.h"
#include "Atmosphere.h"
#include <cmath>

namespace {
    const double g = 9.80665;
    
    struct State {
        double t;
        // Цель
        double Theta_c, v_c, x_gc, y_gc;
        // Носитель
        double Theta_n, v_n, x_gn, y_gn;
        // Ракета
        double Theta_p, v_p, x_gp, y_gp;
    };
    
    // Перегрузки для цели (фрегат — прямолинейное равномерное движение по воде)
    void getTargetOverloads(const State& s, double& n_xa_c, double& n_ya_c) {
        n_xa_c = 0.0;  // постоянная скорость
        n_ya_c = 0.0;  // нет манёвра, движется по прямой y=0
    }
    
    // Перегрузки для носителя (не нужны после отстыковки)
    void getCarrierOverloads(const State& s, double& n_xa_n, double& n_ya_n) {
        n_xa_n = 0.0;
        n_ya_n = 0.0;
    }
    
    // Тангенциальная перегрузка ракеты из таблицы n_xa_r_potr
    double getMissileTanOverload(double t, const R_Data& rdata) {
        const auto& table = rdata.n_xa_r_potr;
        
        if (table.empty()) return 0.0;
        if (t <= table[0][0]) return table[0][1];
        if (t >= table.back()[0]) return table.back()[1];
        
        // Линейная интерполяция
        for (size_t i = 0; i < table.size() - 1; ++i) {
            double t0 = table[i][0];
            double t1 = table[i+1][0];
            double n0 = table[i][1];
            double n1 = table[i+1][1];
            
            if (t >= t0 && t <= t1) {
                if (fabs(t1 - t0) < 1e-12) return n0;
                return n0 + (n1 - n0) * (t - t0) / (t1 - t0);
            }
        }
        return table.back()[1];
    }
    
    // === ЭВТ: пропорциональное сближение с K_g > 1 ===
    double compute_n_ya_p_EWT(const State& s, const GuidanceMethod_Data& m,
                              double epsilon_dot_pc) {
        double k = m.k;           // коэффициент пропорциональности
        double K_g = m.K_g;       // коэффициент компенсации веса (для ЭВТ)
        
        // Формула из пособия (2.26) с K_g вместо 1:
        // n_ya_p^потр = (v_p/g) * k * epsilon_dot + K_g * cos(Theta_p)
        return (s.v_p / g) * k * epsilon_dot_pc + K_g * cos(s.Theta_p);
    }
    
    // === ГОРКА: ЭВТ + манёвр подъёма на конечном участке ===
    double compute_n_ya_p_GORKA(const State& s, const GuidanceMethod_Data& m,
                                 double epsilon_dot_pc, double r_pc,
                                 double r_gorka_start, double n_ya_gorka) {
        if (r_pc > r_gorka_start) {
            // До начала «горки» — обычная ЭВТ
            return compute_n_ya_p_EWT(s, m, epsilon_dot_pc);
        } else {
            // Участок «горки» — дополнительная нормальная перегрузка для подъёма
            // Ракета должна подняться, затем снизиться на цель
            double n_ya_base = compute_n_ya_p_EWT(s, m, epsilon_dot_pc);
            // Добавляем перегрузку подъёма (зависит от текущей фазы горки)
            // Фаза 1: подъём (добавляем положительную перегрузку)
            // Фаза 2: снижение (добавляем отрицательную)
            // Упрощённо: резкий подъём при входе в зону горки
            return n_ya_base + n_ya_gorka * sin(M_PI * (r_gorka_start - r_pc) / r_gorka_start);
        }
    }
}

std::string CalcTrajectory_Euler(
    const CalcParameters& calcParams,
    const InitConditions& inits,
    const LimitConditions& limits,
    const C_Data& cdata,
    const N_Data& ndata,
    const R_Data& rdata,
    const GuidanceMethod_Data& methodData,
    TrajectoryMode mode,
    std::vector<TrajectoryParameters>& result)
{
    double dt = calcParams.dt;
    double r_por = calcParams.r_por;
    double t_max = limits.t_r_max;
    double n_ya_r_max = limits.n_ya_r_max;
    double y_g_r_min = limits.y_g_r_min;
    
    // Параметры для «горки» (извлекаем из methodData или задаём константами)
    double r_gorka_start = methodData.r_gorka_start;  // расстояние до цели для начала горки
    double n_ya_gorka = methodData.n_ya_gorka;        // амплитуда перегрузки горки
    
    // === Инициализация ===
    State s;
    s.t = 0.0;
    
    // Цель
    s.Theta_c = inits.Theta_c0;
    s.v_c = inits.v_c0;
    s.x_gc = inits.x_g_c0;
    s.y_gc = inits.y_g_c0;
    
    // Носитель (начальные = пусковые, но после отстыковки не интегрируем)
    s.Theta_n = inits.Theta_n0;
    s.v_n = inits.v_n0;
    s.x_gn = inits.x_g_pusk;
    s.y_gn = inits.y_g_pusk;
    
    // Ракета
    s.Theta_p = inits.Theta_r0;
    s.v_p = inits.v_r0;
    s.x_gp = inits.x_g_pusk;
    s.y_gp = inits.y_g_pusk;
    
    // Если начальный угол ракеты не задан — направляем в МТВ
    if (s.Theta_p == 0.0 && methodData.auto_theta_r0) {
        double epsilon_pc0 = epsilon(s.x_gp, s.y_gp, s.x_gc, s.y_gc);
        double phi_upr = -asin((s.v_c / s.v_p) * sin(epsilon_pc0 - s.Theta_c));
        s.Theta_p = epsilon_pc0 + phi_upr;
    }
    
    result.clear();
    
    // === Цикл интегрирования ===
    while (s.t <= t_max) {
        // Сохраняем состояние
        TrajectoryParameters tp;
        tp.t = s.t;
        tp.Theta_c = s.Theta_c; tp.v_c = s.v_c; tp.x_g_c = s.x_gc; tp.y_g_c = s.y_gc;
        tp.Theta_n = s.Theta_n; tp.v_n = s.v_n; tp.x_g_n = s.x_gn; tp.y_g_n = s.y_gn;
        tp.Theta_r = s.Theta_p; tp.v_r = s.v_p; tp.x_g_r = s.x_gp; tp.y_g_r = s.y_gp;
        result.push_back(tp);
        
        // === Кинематика ракета-цель ===
        double r_pc = r(s.x_gp, s.y_gp, s.x_gc, s.y_gc);
        double epsilon_pc = epsilon(s.x_gp, s.y_gp, s.x_gc, s.y_gc);
        double r_dot_pc = dot_r(s.x_gp, s.y_gp, s.x_gc, s.y_gc, 
                                     s.v_p, s.v_c, s.Theta_p, s.Theta_c);
        double epsilon_dot_pc = dot_epsilon(s.x_gp, s.y_gp, s.x_gc, s.y_gc,
                                                s.v_p, s.v_c, s.Theta_p, s.Theta_c);
        
        // === Проверка встречи ===
        if (r_pc <= r_por) {
            return "0"; // Успех
        }
        
        // === Перегрузки ===
        double n_xa_c, n_ya_c;
        getTargetOverloads(s, n_xa_c, n_ya_c);
        
        double n_xa_n, n_ya_n;
        getCarrierOverloads(s, n_xa_n, n_ya_n);
        
        double n_xa_p = getMissileTanOverload(s.t, rdata);
        
        double n_ya_p;
        switch (mode) {
            case TrajectoryMode::EWT:
                n_ya_p = compute_n_ya_p_EWT(s, methodData, epsilon_dot_pc);
                break;
            case TrajectoryMode::GORKA:
                n_ya_p = compute_n_ya_p_GORKA(s, methodData, epsilon_dot_pc, 
                                               r_pc, r_gorka_start, n_ya_gorka);
                break;
        }
        
        // === Проверка ограничений ===
        if (fabs(n_ya_p) > n_ya_r_max) {
            return "n_ya_r_max";
        }
        if (s.y_gp < y_g_r_min) {
            return "y_g_r_min";
        }
        
        // === Производные (система ОДУ) ===
        // Цель
        double dTheta_c = (s.v_c != 0) ? ((n_ya_c - cos(s.Theta_c)) * g / s.v_c) : 0;
        double dv_c = (n_xa_c - sin(s.Theta_c)) * g;
        double dx_c = s.v_c * cos(s.Theta_c);
        double dy_c = s.v_c * sin(s.Theta_c);
        
        // Носитель (не интегрируем — отстыкован)
        // double dTheta_n = ...; // пропускаем
        
        // Ракета
        double dTheta_p = (s.v_p != 0) ? ((n_ya_p - cos(s.Theta_p)) * g / s.v_p) : 0;
        double dv_p = (n_xa_p - sin(s.Theta_p)) * g;
        double dx_p = s.v_p * cos(s.Theta_p);
        double dy_p = s.v_p * sin(s.Theta_p);
        
        // === Шаг Эйлера ===
        s.Theta_c += dTheta_c * dt;
        s.v_c += dv_c * dt;
        s.x_gc += dx_c * dt;
        s.y_gc += dy_c * dt;
        
        s.Theta_p += dTheta_p * dt;
        s.v_p += dv_p * dt;
        s.x_gp += dx_p * dt;
        s.y_gp += dy_p * dt;
        
        s.t += dt;
    }
    
    return "t_r_max";
}