set term pdfcairo size 8cm,10cm font "Times,12"
set output 'dof.icub.pdf'

load 'magma.pal'
set key top right
set ylabel 'Density'
set samples 1000

set multiplot layout 4,1
set ytics format "%6.3f"
plot 'mat/DATA_feasible_icub_v2_7_warm_2023-04-30_04-21-27.txt' u ($5):(1/1000.) smooth kdensity ls 1 lw 1 t '7 DOF'
plot 'mat/DATA_feasible_icub_v2_8_warm_2023-04-30_04-21-27.txt' u ($5):(1/1000.) smooth kdensity ls 1 lw 1 t '8 DOF'
plot 'mat/DATA_feasible_icub_v2_9_warm_2023-04-30_04-21-29.txt' u ($5):(1/1000.) smooth kdensity ls 1 lw 1 t '9 DOF'
set xlabel 'Solve time (s)'
plot 'mat/DATA_feasible_icub_v2_10_warm_2023-04-30_04-21-28.txt' u ($5):(1/1000.) smooth kdensity ls 1 lw 1 t '10 DOF'
unset multiplot