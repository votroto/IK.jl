set term pdfcairo size 8cm,8cm font "Times,12"
set output 'infeas.pdf'

load 'magma.pal'
set key top right
set ylabel 'Density'
set samples 1000
set multiplot layout 2,1

# dat/DATA_feas_infeas_icub_v2_7_2023-04-24_04-00-50.txt
# dat/DATA_feas_infeas_kuka_2023-04-24_04-00-51.txt
set title "KUKA LBR iiwa"
#set ytics 1
set xrange [0:3]

plot \
    '<(grep OPTIMAL dat/DATA_uniform_kuka_warm_2023-04-25_19-58-16.txt)' u ($5):(1/679.) smooth kdensity ls 1 lw 2 t 'Feasible', \
    '<(grep OPTIMAL mat/DATA_feasible_kuka_warm_2023-04-30_04-21-27.txt)'         u ($5):(1/1000.) smooth kdensity ls 1 lw 1.5 dt '.' notitle, \
    '<(grep INFEASIBLE dat/DATA_uniform_kuka_warm_2023-04-25_19-58-16.txt)' u ($5):(1/321.) smooth kdensity ls 2 lw 2 t 'Infeasible',\
    '<(grep INFEASIBLE mat/DATA_uniform_kuka_warm_2023-04-30_04-21-28.txt)'       u ($5):(1/335.) smooth kdensity ls 2 lw 1.5 dt '.' notitle, \

set xlabel 'Solve time (s)'
#set ytics 10
set title "iCub"
set xrange [0:1]

plot \
    '<(grep OPTIMAL dat/DATA_feasible_icub_v2_7_warm_2023-04-25_19-58-26.txt)' u ($5):(1/1000.) smooth kdensity ls 1 lw 2 t 'Feasible', \
    '<(grep OPTIMAL mat/DATA_feasible_icub_v2_7_warm_2023-04-30_04-21-27.txt)'    u ($5):(1/1000.) smooth kdensity ls 1 lw 1.5 dt '.' notitle, \
    '<(grep INFEASIBLE dat/DATA_uniform_icub_v2_7_warm_2023-04-25_19-58-16.txt)' u ($5):(1/980.) smooth kdensity ls 2 lw 2 t 'Infeasible',\
    '<(grep INFEASIBLE mat/DATA_uniform_icub_v2_7_warm_2023-04-30_04-21-28.txt)'  u ($5):(1/968.) smooth kdensity ls 2 lw 1.5 dt '.' notitle

unset multiplot



