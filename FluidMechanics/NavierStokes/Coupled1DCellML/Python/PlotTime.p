# Gnuplot script file for plotting data in file "Results/time_x.dat"
# This file is called: plot.p

set terminal postscript portrait enhanced color
set output 'Results/Plot.ps'

set autoscale                          # scale axes automatically
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically
set grid ytics lt 0 lw 3 lc rgb "#bbbbbb"
set grid xtics lt 0 lw 3 lc rgb "#bbbbbb"

set style line 1 lt 1 lw 2 lc rgb "black"
set style line 2 lt 2 lw 2 lc rgb "black"
set style line 3 lt 1 lw 2 lc rgb "red"
set style line 4 lt 2 lw 2 lc rgb "red"
set style line 5 lt 1 lw 2 lc rgb "blue"
set style line 6 lt 2 lw 2 lc rgb "blue"
set style line 7 lt 1 lw 2 lc rgb "green"
set style line 8 lt 2 lw 2 lc rgb "green"
set style line 9 lt 4 lw 2 lc rgb "black"

set multiplot layout 3, 1 title "1DTransient Example"
set key on
#set title "Flow-Distance"
set ylabel "Flow (ml/s)"
set xlabel "Distance (mm)"
#set yrange [-0.2:1.0]
plot 'Results/time_1'  using 1:2 title 'inlet' with lines ls 1,\
     'Results/time_2'  using 1:2 title 'outlet-1' with lines ls 3,\
     'Results/time_3'  using 1:2 title 'outlet-2' with lines ls 5,\
     'Results/time_4'  using 1:2 title 'femoral' with lines ls 7,\
     'Results/time_5'  using 1:2 title 'abd aorta*' with lines ls 2,\
     'Results/time_6'  using 1:2 title 'com iliac*' with lines ls 4,\
     'Results/time_7'  using 1:2 title 'ext iliac*' with lines ls 6,\
     'Results/time_8'  using 1:2 title 'femoral*' with lines ls 8

#set title "Pressure-Distance"
set ylabel "Pressure (mmHg)"
#set yrange [-0.2:1.0]
#set yrange [-0.2:1.0]
plot 'Results/time_1'  using 1:3 title 'inlet' with lines ls 1,\
     'Results/time_2'  using 1:3 title 'outlet-1' with lines ls 3,\
     'Results/time_3'  using 1:3 title 'outlet-2' with lines ls 5,\
     'Results/time_4'  using 1:3 title 'ext iliac' with lines ls 7,\
     'Results/time_5'  using 1:3 title 'femoral' with lines ls 2,\
     'Results/time_6'  using 1:3 title 'deep fem' with lines ls 4,\
     'Results/time_7'  using 1:3 title 'ant tibial' with lines ls 6,\
     'Results/time_8'  using 1:3 title 'pos tibial' with lines ls 8
     
#set title "Conc-Distance"
set ylabel "Conc (mMol)"
set yrange [0.0:1.0]
#set yrange [-0.2:1.0]
plot 'Results/time_1'  using 1:4 title 'Normal' with lines ls 1,\
     'Results/time_2'  using 1:4 title 'Stenosis' with lines ls 3,\
     'Results/time_3'  using 1:4 title 'int iliac' with lines ls 5,\
     'Results/time_4'  using 1:4 title 'ext iliac' with lines ls 7,\
     'Results/time_5'  using 1:4 title 'femoral' with lines ls 2,\
     'Results/time_6'  using 1:4 title 'deep fem' with lines ls 4,\
     'Results/time_7'  using 1:4 title 'ant tibial' with lines ls 6,\
     'Results/time_8'  using 1:4 title 'pos tibial' with lines ls 8
     
unset multiplot
