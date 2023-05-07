set term pdfcairo size 8cm,5cm font "Times,12"
set output 'compare.kuka.pdf'

load 'magma.pal'
set key top right
set xlabel 'Solve time (s)'
set ylabel 'Density'
set samples 1000


set arrow 1 from 4.7,graph(0,0) to 4.7,graph(1,1) nohead lc "black" dt "."
set arrow 2 from 0.44,graph(0,0) to 0.44,graph(1,1) nohead lc "black" dt "."

plot \
    	'dat/METHOD_ltr_2023-04-27_03-06-48.txt' u ($5):(1/100.) smooth kdensity ls 1 lw 1.5 t 'LR', \
	'dat/METHOD_mat_2023-04-27_03-06-47.txt' u ($5):(1/100.) smooth kdensity ls 2 lw 1.5 t 'Mat.', \
	'dat/METHOD_tree_2023-04-27_03-06-47.txt' u ($5):(1/100.) smooth kdensity ls 3 lw 1.5 t 'Par.'