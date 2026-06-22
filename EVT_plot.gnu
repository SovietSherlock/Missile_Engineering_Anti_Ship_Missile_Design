# Gnuplot script for ЭВТ ракеты
# Usage: gnuplot -persist EVT_plot.gnu

set terminal wxt size 800,600 enhanced font 'Arial,12'
set title 'ЭВТ ракеты - Траектория полёта'
set xlabel 'x_g, м'
set ylabel 'y_g, м'
set grid
plot 'EVT_trajectory.gra' using 1:2 with lines lw 2 lt rgb 'red' title 'Ракета (ЭВТ)'
pause -1 'Нажмите Enter для следующего графика...'

set title 'ЭВТ ракеты - Скорость ракеты'
set xlabel 't, с'
set ylabel 'v, м/с'
set grid
plot 'EVT_velocity.gra' using 1:2 with lines lw 2 lt rgb 'blue' title 'v(t)'
pause -1 'Нажмите Enter для следующего графика...'

set title 'ЭВТ ракеты - Координаты'
set xlabel 't, с'
set ylabel 'Координата, м'
set grid
plot 'EVT_x_t.gra' using 1:2 with lines lw 2 lt rgb 'green' title 'x_g(t)', \
     'EVT_y_t.gra' using 1:2 with lines lw 2 lt rgb 'red' title 'y_g(t)'
pause -1 'Нажмите Enter для следующего графика...'

set title 'ЭВТ ракеты - Угол траектории'
set xlabel 't, с'
set ylabel 'Theta, град'
set grid
plot 'EVT_theta.gra' using 1:3 with lines lw 2 lt rgb 'purple' title 'Theta(t)'
pause -1 'Нажмите Enter для следующего графика...'

set title 'ЭВТ ракеты - Нормальная перегрузка'
set xlabel 't, с'
set ylabel 'n_ya'
set grid
plot 'EVT_n_ya.gra' using 1:2 with lines lw 2 lt rgb 'orange' title 'n_ya(t)'
pause -1 'Нажмите Enter для завершения...'
