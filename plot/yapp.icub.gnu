set term pdfcairo size 7cm,5cm font "Times,12"
set output 'cloud.globloc.icub8.pdf'

ntics = 3.0

function $rround(x) << EOF
      return ceil(x*10)/10.0;
EOF

stats 'CLOUD_feasible_icub_v2_8_2024-02-23_19-59-56.txt' using 1 name 'x' nooutput
stats 'CLOUD_feasible_icub_v2_8_2024-02-23_19-59-56.txt' using 2 name 'y' nooutput
stats 'CLOUD_feasible_icub_v2_8_2024-02-23_19-59-56.txt' using 3 name 'z' nooutput

set xtics $rround((($rround(x_max) - $rround(x_min))/ntics))
set ytics $rround((($rround(y_max) - $rround(y_min))/ntics))
set ztics $rround((($rround(z_max) - $rround(z_min))/ntics))

set xrange [$rround(x_min):$rround(x_max)]
set yrange [$rround(y_min):$rround(y_max)]
set zrange [$rround(z_min):$rround(z_max)]

set xlabel "x"
set ylabel "y"
set zlabel "z"

set palette defined (-0.6 '#b2182b', -0.1 '#b2182b', 0 '#fddbc7', 0.1 '#2166ac', 0.3 '#2166ac')

set cbrange [-0.6:0.3]
set cbtics 0.4

set xyplane 0
set view 60, 160, 1, 1.2
splot "CLOUD_feasible_icub_v2_8_2024-02-23_19-59-56.txt" u 1:2:3:($5 - $8) w points palette pt 7 ps 0.1 notitle