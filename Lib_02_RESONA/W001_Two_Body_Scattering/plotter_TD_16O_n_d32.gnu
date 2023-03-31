
set rmargin 0.3
set tmargin 0.3
set size 0.8, 0.8

set xlabel "r  [fm]"
set ylabel "[MeV]"
p[0:23][-18:12] \
"fort.2001" u 1:2  w l ti "Total V(r)", \
"fort.2001" u 1:5  w l ti "Woods-Saxon", \
"fort.2001" u 1:3  w l ti "Centrifugal", \
"fort.2011" u 1:3      ti "Conf. Pot.", \
"fort.2011" u 1:(($5)*20)  w l lc "black" dt 3 ti "20*{/Symbol r}_{n}(t) at ct = 0"
pause  1

set terminal postscript eps enhanced color font "Century, 20"
#set output "hogee.eps"
set output "| epstopdf -f -o=hogee.pdf"
replot
unset output
set terminal x11
unset label
reset


set xlabel "ct  [fm]"
set ylabel "[1]"
p[][] "fort.2021" u 1:3 w l ti "Gamma(t)"
pause  1
p[][] "fort.2021" u 1:2 w l ti "P_{surv}(t)"
pause  1


set xlabel "r  [fm]"
p[0:40][0:0.32] "fort.2500" u 1:2 ti "{/Symbol r}_{n}(t) at ct = 0", \
"fort.2011" u 1:5  w l    ti ""
pause  1
p[0:40][0:0.32] "fort.2501" u 1:2 ti "ct = 2000 fm"
pause 0.4
p[0:40][0:0.32] "fort.2502" u 1:2 ti "ct = 4000 fm"
pause 0.4
p[0:40][0:0.32] "fort.2503" u 1:2 ti "ct = 6000 fm"


