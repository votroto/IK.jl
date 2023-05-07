set term pdfcairo size 8cm,5cm font "Times,12"
set output 'times.kuka.pdf'

load 'magma.pal'
set key top right
set xrange [0 to 15]
set xlabel 'Solve time (s)'
set ylabel 'Density'
set samples 1000


set arrow 1 from 4.7,graph(0,0) to 4.7,graph(1,1) nohead lc "black" dt "."
set arrow 2 from 0.44,graph(0,0) to 0.44,graph(1,1) nohead lc "black" dt "."

plot \
    'dat/DATA_feasible_kuka_warm_2023-04-25_19-58-16.txt' u ($5):(1/1000.) smooth kdensity ls 1 lw 1.5 t '(Gurobi) QCQP', \
    'dat/DATA_SCIP_bench_kuka_2023-04-24_03-10-49.txt' u ($5):(1/1000.) smooth kdensity t '(SCIP) QCQP' ls 1 lw 1.5 dt "-", \
    'dat/tims.txt' u ($1+$2):(1/300.) smooth kdensity t 'SOS' ls 2 lw 1.5 \
