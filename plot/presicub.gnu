set term pdfcairo size 8cm,6cm font "Times,12"
set output 'pres.icub.pdf'

set palette viridis
set key top right
set ylabel 'Density'
set samples 1000

set multiplot layout 2,1
set ytics 1
set title "iCub 7-DOF"
plot \
    'cluster/DATA_feasible_icub_v2_7_angdiff_2024-04-23_11-26-03.txt' u ($5):(1/200.) smooth cnorm ls 1 lw 1 t 'angdiff', \
    'cluster/DATA_feasible_icub_v2_7_feas_2024-04-23_11-26-03.txt' u ($5):(1/200.) smooth cnorm ls 3 lw 1 t 'feasible', \
    'cluster/DATA_feasible_icub_v2_7_l2_2024-04-23_11-26-04.txt' u ($5):(1/200.) smooth cnorm ls 7 lw 1 t 'squared ℓ_2' 
set title "iCub 10-DOF"
set xlabel 'Solve time (s)'
#plot \
#    'cluster/DATA_feasible_icub_v2_10_angdiff_2024-04-23_11-26-03.txt' u ($5):(1/200.) smooth cnorm ls 1 lw 1 t 'angdiff', \
#    'cluster/DATA_feasible_icub_v2_10_feas_2024-04-23_11-59-56.txt' u ($5):(1/200.) smooth cnorm ls 3 lw 1 t 'feasible', \
#    'cluster/DATA_feasible_icub_v2_10_l2_2024-04-23_11-33-12.txt' u ($5):(1/200.) smooth cnorm ls 7 lw 1 t 'squared ℓ_2'

plot \
    'local/DATA_feasible_icub_v2_10_angdiff_2024-04-23_12-02-46.txt' u ($5):(1/200.) smooth cnorm ls 1 lw 1 t 'angdiff', \
    'local/DATA_feasible_icub_v2_10_feas_2024-04-23_12-02-56.txt' u ($5):(1/200.) smooth cnorm ls 3 lw 1 t 'feasible', \
    'local/DATA_feasible_icub_v2_10_l2_2024-04-23_12-02-30.txt' u ($5):(1/200.) smooth cnorm ls 7 lw 1 t 'squared ℓ_2'

unset multiplot