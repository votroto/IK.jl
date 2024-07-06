set term pdfcairo size 8cm,5cm font "Times,12"
set output 'times.postgb.pdf'

load 'plot/magma.pal'
set key bottom right
set xlabel 'Solve time (s)'
set ylabel 'CDF'
set xrange [0:20]
plot \
    'DATA_gurobi.txt' u ($5):(1) smooth cnorm ls 1 lw 1.5 t 'Gurobi', \
    'DATA_SCIP2.txt' u ($5):(1) smooth cnorm t 'SCIP' ls 1 lw 1.5 dt "-", \
    'execution_times_scaled.txt' u ($2+(1-$1)*100):(1) smooth cnorm t 'SOS' ls 2 lw 1.5 \
