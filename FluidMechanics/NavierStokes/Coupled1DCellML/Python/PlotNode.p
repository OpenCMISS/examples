# Gnuplot script file for plotting data in file "Results/node_x.dat"
# This file is called: plot.p
set termopt enhanced
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
set style line 9 lt 1 lw 2 lc rgb "orange"
set style line 10 lt 2 lw 2 lc rgb "orange"
set style line 11 lt 1 lw 2 lc rgb "purple"
set style line 12 lt 2 lw 2 lc rgb "purple"
set style line 13 lt 1 lw 2 lc rgb "cyan"
set style line 14 lt 2 lw 2 lc rgb "cyan"
set style line 15 lt 1 lw 2 lc rgb "yellow"
set style line 16 lt 2 lw 2 lc rgb "yellow"
set style line 17 lt 1 lw 2 lc rgb "navy"
set style line 18 lt 2 lw 2 lc rgb "navy"

#set logscale y
#set y2tics
#set ytics nomirror
#set y2label "Conc (mMol)"
#set y2range [0:1]

set multiplot layout 3, 1 title "1DTransient Example"
set key on
#set key samplen 2 spacing .5 font ",8"
#set title "Flow-Time"
set ylabel "Flow (ml.s^{-1})"
set xlabel "Time (s)"
#set yrange [0.0001:0.001]
plot 'Results/node_1'  using 1:2 title 'Normal' with lines ls 1,\
     'Results/node_2'  using 1:2 title 'A = 2x' with lines ls 3,\
     'Results/node_3'  using 1:2 title 'A = 3x' with lines ls 5,\
     'Results/node_4'  using 1:2 title 'A = 4x' with lines ls 7,\
     'Results/node_5'  using 1:2 title 'A = 5x' with lines ls 9,\
     'Results/node_6'  using 1:2 title 'A = 6x' with lines ls 11,\
     'Results/node_7'  using 1:2 title 'A = 7x' with lines ls 13,\
     'Results/node_8'  using 1:2 title 'A = 8x' with lines ls 15,\
     'Results/node_9'  using 1:2 title 'A = 9x' with lines ls 17,\
     'Results/node_10' using 1:2 title 'A = 10x' with lines ls 10

#set title "Pressure-Time"
set ylabel "Pressure (mmHg)"
#set yrange [-0.2:1.0]
#set yrange [-0.2:1.0]
plot 'Results/node_1'  using 1:3 title 'Normal' with lines ls 1,\
     'Results/node_2'  using 1:3 title 'A = 2x' with lines ls 3,\
     'Results/node_3'  using 1:3 title 'A = 3x' with lines ls 5,\
     'Results/node_4'  using 1:3 title 'A = 4x' with lines ls 7,\
     'Results/node_5'  using 1:3 title 'A = 5x' with lines ls 9,\
     'Results/node_6'  using 1:3 title 'A = 6x' with lines ls 11,\
     'Results/node_7'  using 1:3 title 'A = 7x' with lines ls 13,\
     'Results/node_8'  using 1:3 title 'A = 8x' with lines ls 15,\
     'Results/node_9'  using 1:3 title 'A = 9x' with lines ls 17,\
     'Results/node_10' using 1:3 title 'A = 10x' with lines ls 10
     
#set title "Conc-Time"
set ylabel "Conc (mMol)"
#set yrange [-0.2:1.0]
#set yrange [-0.2:1.0]
plot 'Results/node_1'  using 1:4 title 'Arteries' with lines ls 1,\
     'Results/node_2'  using 1:4 title 'Veins' with lines ls 3,\
     'Results/node_3'  using 1:4 title 'Kidney' with lines ls 5,\
     'Results/node_4'  using 1:4 title 'ext iliac' with lines ls 7,\
     'Results/node_5'  using 1:4 title 'femoral' with lines ls 9,\
     'Results/node_6'  using 1:4 title 'deep fem' with lines ls 2,\
     'Results/node_7'  using 1:4 title 'ant tibial' with lines ls 4,\
     'Results/node_8'  using 1:4 title 'pos tibial' with lines ls 6,\
     'Results/node_9'  using 1:4 title 'pos tibial' with lines ls 8,\
     'Results/node_10' using 1:4 title 'pos tibial' with lines ls 10
     
unset multiplot
