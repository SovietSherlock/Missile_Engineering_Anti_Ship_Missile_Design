#include "GuidanceMethods.h"
#include "KinematicFunc2D.h"

#include <cmath>

#include "LibConstFunc.h"

/// Функция для вычисления потребной перегрузки по методу чистой погони (самонаведение)
double n_ya_potr_PurePursuit(double dot_epsilon_rc, double v_r, double Theta_r)
{
    return v_r / g * dot_epsilon_rc + cos(Theta_r);
}

/// Функция для вычисления потребной перегрузки по методу наведения по лучу (теленаведение), для общего случая
double n_ya_potr_BeamRiding(double epsilon_l, double dot_epsilon_l, double ddot_epsilon_l, double r_nr, double dot_r_nr, double Theta_r, double Theta_n, double a_xa_r, double a_xa_n, double a_ya_n)
{
    return 1.0 / (g * cos(epsilon_l - Theta_r)) * (ddot_epsilon_l * r_nr + 2.0 * dot_epsilon_l * dot_r_nr +
                                                   a_xa_n * sin(Theta_n - epsilon_l) + a_ya_n * cos(epsilon_l - Theta_n) -
                                                   a_xa_r * sin(Theta_r - epsilon_l)) + cos(Theta_r);
}

double n_ya_potr_Proportional(double k, double dot_epsilon_rc, double v_r, double Theta_r)
{
    return v_r/ g * k * dot_epsilon_rc + cos(Theta_r);
}

double n_ya_potr_EVT(double v_r, double dot_epsilon_rc, double Theta_r, double K_g, double coeff)
{
    double result = (v_r / g) * dot_epsilon_rc * coeff + K_g * cos(Theta_r);
    return result;
}

/// Энергетически выгодная траектория с наведением на виртуальную цель
double n_ya_potr_EVT_virtual(
    double v_r, double v_c,
    double Theta_r, double Theta_c,
    double dot_epsilon_virtual,
    double K_g, double k)
{
    // Классическая ЭВТ, но с угловой скоростью на виртуальную цель
    return (v_r / g) * k * dot_epsilon_virtual + K_g * cos(Theta_r);
}

/// Маршевый полёт (закон управления высотой)
double n_ya_potr_March(
    double v_r, double Theta_r,
    double y_g_r,
    double H_march,
    double K_H_march,
    double K_v_march)
{
    // n_ya = K_H*(H_марш - y) - K_v*v*sin(Theta) + cos(Theta)
    return K_H_march * (H_march - y_g_r)
         - K_v_march * v_r * sin(Theta_r)
         + cos(Theta_r);
}

/// Подъём на высоту "горки"
double n_ya_potr_GorkaClimb(
    double v_r, double Theta_r,
    double y_g_r,
    double H_gorka,
    double K_H_gorka,
    double K_v_gorka)
{
    // n_ya = K_H_горки*(H_горки - y) - K_v_горки*v*sin(Theta) + cos(Theta)
    return K_H_gorka * (H_gorka - y_g_r)
         - K_v_gorka * v_r * sin(Theta_r)
         + cos(Theta_r);
}

/// Главная функция управления многофазной траекторией
double n_ya_potr_TwoPhase(
    double v_r, double Theta_r,
    double y_g_r,
    double r_rc,
    double epsilon_rc,
    double v_c, double Theta_c,
    double dot_epsilon_virtual,
    double dot_epsilon_rc,
    const GuidanceMethod_Data& mtd,
    int& phase,
    double n_ya_max)
{
    // ========================================================================
    // ФАЗА 0: ЭНЕРГЕТИЧЕСКИ ВЫГОДНАЯ ТРАЕКТОРИЯ (на виртуальную цель)
    // ========================================================================
    if (phase == 0)
    {
        // Перегрузка по ЭВТ на виртуальную цель
        double n_ya_evt = n_ya_potr_EVT_virtual(
            v_r, v_c,
            Theta_r, Theta_c,
            dot_epsilon_virtual,
            mtd.K_g,
            mtd.k
        );
        
        // Перегрузка для маршевого режима (для сравнения)
        double n_ya_march = n_ya_potr_March(
            v_r, Theta_r,
            y_g_r,
            mtd.H_march,
            mtd.K_H_march,
            mtd.K_v_march
        );
        
        // УСЛОВИЕ ПЕРЕХОДА: 
        // ЭВТ активна при: Θ_p >= 0
        // Переход на марш происходит когда:
        // 1) Θ_p < 0 (ракета перешла в пикирование)
        // 2) -n_ya_max <= n_ya_march <= n_ya_max (маршевая перегрузка в допустимом диапазоне)
        // ================================================================
        if (Theta_r < 0.0 && n_ya_march < n_ya_max && n_ya_march > -n_ya_max)
        {
            phase = 1;  // переходим на маршевый участок
            return n_ya_march;
        }
        return n_ya_evt;
    }
    
    // ========================================================================
    // ФАЗА 1: МАРШЕВЫЙ ПОЛЁТ
    // ========================================================================
    if (phase == 1)
    {
        double n_ya = n_ya_potr_March(
            v_r, Theta_r,
            y_g_r,
            mtd.H_march,
            mtd.K_H_march,
            mtd.K_v_march
        );
        
        
        // УСЛОВИЕ ПЕРЕХОДА: достигли безопасного расстояния для горки
        if (r_rc <= mtd.r_save)
        {
            phase = 2;  // переходим к подъёму
        }
        
        return n_ya;
    }
    
    // ========================================================================
    // ФАЗА 2: ПОДЪЁМ НА ВЫСОТУ "ГОРКИ" И ПОЛЕТ НА НЕЙ
    // ========================================================================
    if (phase == 2)
    {
        double r_crit_calc = calc_r_critical(
            v_r, v_c,
            mtd.H_gorka,
            mtd.k_dive,
            n_ya_max
        );

        // УСЛОВИЕ ПЕРЕХОДА: достигли критического расстояния для пикирования
        if (r_rc <= r_crit_calc && r_crit_calc > 0.0|| r_rc < 750.0)
        {
            phase = 3;  
        }
        
        double n_ya = n_ya_potr_GorkaClimb(
            v_r, Theta_r,
            y_g_r,
            mtd.H_gorka,
            mtd.K_H_gorka,
            mtd.K_v_gorka
        );
        //Условие непревышения максималльной перегрузки на подъеме на высоту горки
        if (n_ya >= n_ya_max)      return n_ya_max;
        if (n_ya <= -n_ya_max)     return -n_ya_max;
        return n_ya;
    }

    // ========================================================================
    // ФАЗА 3: ПИКИРОВАНИЕ (для визуализации/логирования)
    // ========================================================================
    if (phase == 3)
    {
        return n_ya_potr_Proportional(mtd.k_dive, dot_epsilon_rc, v_r, Theta_r);
    }
}