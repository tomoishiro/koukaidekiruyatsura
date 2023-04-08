set tmargin 0.3
set rmargin 0.3

set ylabel "R_{max}  [fm]"
set xlabel "E_{SP}  [MeV]"

e0 = 0.415071
g0 = 0.0166812
c1 = -7.45974
c2 = -96.022
c3 =  118.657
f(x) = c1*atan( (x-e0)/(0.5*g0) )   +c2*x +c3

fit[0.38:0.45]  f(x)  "Fine.dat" u 4:1  via  e0,g0,c1,c2,c3

#set xtics 2.0, 0.1, 10.0
p [][] "Fine.dat"  u 4:1 w p lt 6 lc "black"  ti "", \
f(x) w l lt 3 lc "red" ti "Fitted"
pause -1

set size 0.6, 0.6
set terminal postscript eps enhanced color font "Century, 20"
set output "| epstopdf -f -o=hogee_stabfit.pdf"
replot
unset output
set terminal x11

