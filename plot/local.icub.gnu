set term pdfcairo size 8cm,5cm font "Times,12"
set output 'local.icub.pdf'

load 'magma.pal'
set key top right
set ylabel 'Density'
set samples 1000
set xrange [-1 to 1]
set xlabel 'Relative gap'

plot \
	'dat/DATA_feasible_icub_v2_7_warm_2023-04-25_19-58-26.txt'  u (($4-$3)/$3):(1/1000.) smooth kdensity ls 1 lw 1 t '7 DOF', \
	'dat/DATA_feasible_icub_v2_8_warm_2023-04-25_19-58-16.txt'  u (($4-$3)/$3):(1/1000.) smooth kdensity ls 2 lw 1 t '8 DOF', \
	'dat/DATA_feasible_icub_v2_9_warm_2023-04-25_19-58-16.txt'  u (($4-$3)/$3):(1/1000.) smooth kdensity ls 3 lw 1 t '9 DOF', \
	'dat/DATA_feasible_icub_v2_10_warm_2023-04-25_19-58-16.txt' u (($4-$3)/$3):(1/303.) smooth kdensity ls 4 lw 1 t '10 DOF'
	#'dat/DATA_feasible_kuka_warm_2023-04-25_19-58-16.txt'  u (($4-$3)/$3):(1/1000.) smooth kdensity ls 5 lw 1 dt '-' t 'KUKA'
unset multiplot