# Gnuplot script file for plotting data in file "Result_nodex.dat"
# This file is called: plot.p

set terminal postscript portrait enhanced color
set output 'Result.ps'

set autoscale                          # scale axes automatically
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set grid ytics lt 0 lw 3 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 3 lc rgb "#bbbbbb"
#set xrange [0:0.8]
set style line 1 lt 1 lw 1 lc rgb "blue"
set style line 2 lt 2 lw 1 lc rgb "blue"
set style line 3 lt 1 lw 1 lc rgb "red"
set style line 4 lt 2 lw 1 lc rgb "red"
set style line 5 lt 1 lw 1 lc rgb "green"
set style line 6 lt 2 lw 1 lc rgb "green"
set style line 7 lt 1 lw 1 lc rgb "black"
set style line 8 lt 2 lw 1 lc rgb "black"

set multiplot layout 3, 1 title "1DTransient Example"
set key off
#set title "Flow-Time"
set ylabel "Flow (ml/s)"
#set yrange [-0.2:1.0]
plot "Result_node5"  using 1:2 title 'dsc Aorta' with lines ls 1,\
     "Result_node9"  using 1:2 title 'dsc Aorta' with lines ls 3,\
     "Result_node31"  using 1:2 title 'com Iliac*' with lines ls 5,\
     "Result_node87"  using 1:2 title 'femoral*' with lines ls 7

set key off
#set title "Pressure-Time"
set ylabel "Pressure (mmHg)"
#set yrange [-0.2:1.0]
plot "Result_node5"  using 1:3 title 'dsc Aorta' with lines ls 1,\
     "Result_node9"  using 1:3 title 'dsc Aorta' with lines ls 3,\
     "Result_node31"  using 1:3 title 'com Iliac*' with lines ls 5,\
     "Result_node87"  using 1:3 title 'femoral*' with lines ls 7

set key off
#set title "Area-Time"
set xlabel "Time (s)"
set ylabel "Area (mm2)"
#set yrange [-0.2:1.0]
plot "Result_node5"  using 1:4 title 'dsc Aorta' with lines ls 1,\
     "Result_node9"  using 1:4 title 'dsc Aorta' with lines ls 3,\
     "Result_node31"  using 1:4 title 'com Iliac*' with lines ls 5,\
     "Result_node87"  using 1:4 title 'femoral*' with lines ls 7

unset multiplot
