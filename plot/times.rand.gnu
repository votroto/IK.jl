set terminal pdf enhanced size 8cm,5cm font "Times,12"
set output 'times.generic.pdf'

load 'magma.pal'
set key top right
set ylabel 'Density'
set samples 2000
set xrange [0.1:104]
set xlabel 'Solve time (s)'
set logscale x 

plot \
	'dat/DATA_feasible_rand_orth_warm_2023-04-25_19-58-25.txt' u ($5):(1/1000.) smooth kdensity ls 1 lw 2 t 'Orth.', \
	'dat/DATA_feasible_rand_orth_cold_2023-04-25_19-58-16.txt' u ($5):(1/1000.) smooth kdensity ls 1 lw 1.5 dt '.' t '(cold) Orth.', \
	'dat/DATA_feasible_rand_4rad_warm_2023-04-25_19-58-16.txt' u ($5):(1/1000.) smooth kdensity ls 2 lw 2 t '4 rad', \
	'dat/DATA_feasible_rand_4rad_cold_2023-04-25_19-58-16.txt' u ($5):(1/1000.) smooth kdensity ls 2 lw 1.5 dt '.' t '(cold) 4 rad', \
	'dat/DATA_feasible_rand_6rad_warm_2023-04-25_19-58-40.txt' u ($5):(1/1000.) smooth kdensity ls 5 lw 2 t '6 rad', \
	'dat/DATA_feasible_rand_6rad_cold_2023-04-25_19-58-16.txt' u ($5):(1/1000.) smooth kdensity ls 5 lw 1.5 dt '.' t '(cold) 6 rad', \

#plot \
#	'mat/DATA_feasible_rand_orth_warm_2023-04-30_04-21-27.txt' u ($5):(1/1000.) smooth kdensity ls 1 lw 2 t 'Orth.', \
#	'mat/DATA_feasible_rand_orth_cold_2023-04-30_04-21-27.txt' u ($5):(1/1000.) smooth kdensity ls 1 lw 1.5 dt '.' t '(cold) Orth.', \
#	'mat/DATA_feasible_rand_4rad_warm_2023-04-30_04-21-27.txt' u ($5):(1/1000.) smooth kdensity ls 2 lw 2 t '4 rad', \
#	'mat/DATA_feasible_rand_4rad_cold_2023-04-30_04-21-27.txt' u ($5):(1/1000.) smooth kdensity ls 2 lw 1.5 dt '.' t '(cold) 4 rad', \
#	'mat/DATA_feasible_rand_6rad_warm_2023-04-30_04-21-28.txt' u ($5):(1/1000.) smooth kdensity ls 5 lw 2 t '6 rad', \
#	'mat/DATA_feasible_rand_6rad_cold_2023-04-30_04-21-27.txt' u ($5):(1/1000.) smooth kdensity ls 5 lw 1.5 dt '.' t '(cold) 6 rad', \



# cat dat/DATA_feasible_rand_orth_warm_2023-04-25_19-58-25.txt | datamash -C -t' ' --format '%.2e' mean 1 mean 2 mean 5
