set terminal pdf enhanced size 8cm,5cm font "Times,12"
set output 'times.gb.rand.pdf'

load 'plot/magma.pal'
set key top left
set ylabel 'CDF'
set samples 2000
set xlabel 'Solve time (s)'

plot \
	'ctimes' u ($1+$2):(1) smooth cnorm ls 1 lw 2 t 'Truncated GB', \
