set term pdfcairo size 7cm,5cm font "Times,12"
set output 'compare.kuka.pdf'

set ztics 0.5, 0.5
set xtics 0.4
set ytics 0.4

set zrange [0:1]

set xlabel "x"
set ylabel "y"
set zlabel "z"

load 'magma.pal'
set cbrange [0:2]
unset colorbox
set xyplane 0
set view 60, 160, 1, 1.2
splot "< awk '{if($6 == \"OPTIMAL\") print}' DATA_kuka_yapp_2024-02-21_20-02-01.txt" u 1:2:3:5 w points palette pt 5 ps 0.1 notitle, \
      "< awk '{if($6 == \"INFEASIBLE\") print}' DATA_kuka_yapp_2024-02-21_20-02-01.txt" u 1:2:3 w points pt 5 ps 0.01 lc rgb "#33aaaaaa" notitle