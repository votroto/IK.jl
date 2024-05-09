set term pdfcairo transparent size 8cm,4cm font "Times,12"
set output 'pres.canadarm.pdf'

set palette viridis
set key top right
set ylabel 'Density'
set samples 1000
set ytics 1

set xlabel 'Solve time (s)'

plot \
    'cluster/DATA_feasible_canadarm_angdiff_2024-04-23_11-26-03.txt' u ($5):(1/200.) smooth cnorm ls 1 lw 1 t 'angdiff', \
    'cluster/DATA_feasible_canadarm_feas_2024-04-23_11-26-04.txt' u ($5):(1/200.) smooth cnorm ls 3 lw 1 t 'feasible', \
    'cluster/DATA_feasible_canadarm_l2_2024-04-23_11-26-03.txt' u ($5):(1/200.) smooth cnorm ls 7 lw 1 t 'squared ℓ_2'


unset multiplot