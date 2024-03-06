# csvを読み込んでプロット

# 等圧線
set datafile separator ","
set xrange [-10:30]
set yrange [-10:10]
set table 'colormap.dat'
splot "MAC_example_5000.csv"
unset table

set contour base
set cntrparam level incremental -1.5, 0.05, 1.5
unset surface
set table 'contour.dat'
splot "MAC_example_5000.csv"
unset table

reset 
set xrange [-10:30]
set yrange [-10:10]
set size ratio 0.5
unset key
# set palette defined (0 '#0fffee',1 '#0090ff', 2 '#000fff',3 '#000090',4 '#ffffff',5 '#7f0000', 6 '#ee0000', 7 '#ff7000', 8 '#ffee00')
set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#ffffff',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
plot 'colormap.dat' with image, 'contour.dat' with lines lt -1 lw 1

# 等渦度線
set datafile separator ","
set xrange [-10:30]
set yrange [-10:10]
set table 'colormap.dat'
splot "MAC_example_omega.csv"
unset table

set contour base
set cntrparam level incremental -7.1, 0.2, 7.1
unset surface
set table 'contour.dat'
splot "MAC_example_omega.csv"
unset table

reset 
set xrange [-10:30]
set yrange [-10:10]
set size ratio 0.5
unset key
# set palette defined (0 '#0fffee',1 '#0090ff', 2 '#000fff',3 '#000090',4 '#ffffff',5 '#7f0000', 6 '#ee0000', 7 '#ff7000', 8 '#ffee00')
set palette defined ( 0 '#000090',1 '#000fff',2 '#0090ff',3 '#0fffee',4 '#ffffff',5 '#ffee00',6 '#ff7000',7 '#ee0000',8 '#7f0000')
plot 'colormap.dat' with image, 'contour.dat' with lines lt -1 lw 1

# 無次元係数
set datafile separator ","
set yrange [-2.5:2.5]
set xlabel "time[s]"
set key font "Arial, 12"
set key box
plot "history_dimless_coef.csv" using 1:2 with lines lw 2 title "CD",\
"history_dimless_coef.csv" using 1:3 with lines lw 2 title "CL",\
"history_dimless_coef.csv" using 1:4 with lines lw 2 title "Cp1",\
"history_dimless_coef.csv" using 1:5 with lines lw 2 title "Cp2"

# 誤差変化
set datafile separator ","
set yrange [1e-6:1e-3]
set logscale y
set format y "1e%T"
set xlabel "time[s]"
set ylabel "RMS error for p"
plot "history_of_rhs2.csv" using 1:2 with lines lw 2
unset key
set xlabel "time[s]"
set ylabel "RMS error for p"
replot