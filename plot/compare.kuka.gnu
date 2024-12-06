set term pdfcairo size 8cm,5cm font "Times,12"
set output 'kdensity.kuka.pdf'

set key top right
set xlabel 'Solve time (s)'
set ylabel 'Density'
set samples 1000


plot 'DATA_feasible_kuka_2024-12-04_22-48-58.txt' u ($5):(1/100.) smooth kdensity ls 1 lw 1.5 t 'Feasible Kuka'