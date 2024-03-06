# wave equations
reset
set xrange[0.0 : 10.0]
set yrange[-1.0 : 1.5]
set pointsize 1.5
set xlabel "x" font "Arial, 15"
set ylabel "u" font "Arial, 15"
set xtics 2
set tics font "Arial, 15"
set key top right font "Arial, 15"

set datafile separator ","
plot "u.csv" using 1:2 with lines lw 2 title "numerical",\
"u.csv" using 1:3 with lines lw 2 title "initial",\
"u.csv" using 1:4 with lines lw 2 title "exact"

# burgers equations
reset
set xrange[0.0 : 200.0]
set yrange[-0.5 : 2.5]
set pointsize 1.5
set xlabel "x" font "Arial, 15"
set ylabel "u" font "Arial, 15"
# set xtics 10
set tics font "Arial, 15"
set key top right font "Arial, 15"

set datafile separator ","
plot "u_muscl_minmod_final.csv" using 1:3 with lines lw 2 title "initial",\
"u_muscl_minmod_final.csv" using 1:4 with lines lw 2 title "exact",\
"u_upw1_final.csv" using 1:2 with lines lw 2 title "upwind1",\
"u_muscl_minmod_copy.csv" using 1:2 with lines lw 2 title "minmod"

# see transitions
reset
set xrange[0.0 : 200.0]
set yrange[-0.5 : 2.5]
set pointsize 1.5
set xlabel "x" font "Arial, 15"
set ylabel "u" font "Arial, 15"
# set xtics 10
set tics font "Arial, 15"
set key top right font "Arial, 15"

set datafile separator ","
plot "u_upw1_trans.csv" using 1:2 with lines lw 2 title "upw1:step=100",\
"u_upw1_final.csv" using 1:2 with lines lw 2 title "upw1:step=200",\
"u_muscl_minmod_trans.csv" using 1:2 with lines lw 2 title "minmod:step=100",\
"u_muscl_minmod_final.csv" using 1:2 with lines lw 2 title "minmod:step=200"