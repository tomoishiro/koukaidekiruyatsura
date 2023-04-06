
set rmargin 0.6
set tmargin 0.6
set size 0.8, 0.8
set key spacing 1.6

set xlabel "r  [fm]"
set ylabel "[MeV]"
p[0:23][-18:12] \
"fort.2001" u 1:2  w l ti "Total V(r)", \
"fort.2001" u 1:5  w l ti "Woods-Saxon", \
"fort.2001" u 1:3  w l ti "Centrifugal", \
"fort.2011" u 1:3      ti "Conf. Pot.", \
"fort.2011" u 1:(($5)*20)  w l lc "black" dt 3 ti "20*{/Symbol r}_{n}(t) at ct = 0"
pause  1


set xlabel "ct  [fm]"
set ylabel "{/Symbol G}(t)  [MeV]"
p[0:6200][0:] "fort.2021" u 1:3 w l ti ""
pause  1
set ylabel "[1]"
p[][] "fort.2021" u 1:2 w l ti "P_{surv}(t)"
pause  1
#----------------------------------------------FITTING
gam = 0.017 #MeV
d0 = 0.03
usr(x) = exp(-gam*x/197.032) + d0
   fit[3000:5000] usr(x)  "fort.2021" u  1:2  via gam, d0
set xlabel  "ct  [fm]"
set ylabel  "P_{surv}(t)"
p[0:5400][:] \
"fort.2021" u  1:(($2))   w l  dt 2 lc 8 lw 3  ti "^{16}O+n, P_{surv}(t)", \
            usr(x)        w l  dt 1 lc 1 lw 3  ti "Fitted"
pause 11

set terminal postscript eps enhanced color font "Century, 20"
#set output "hogee.eps"
set output "| epstopdf -f -o=hogee.pdf"
replot
unset output
set terminal x11
unset label
reset


set xlabel "r  [fm]"
p[0:40][0:0.32] "fort.2500" u 1:2 ti "{/Symbol r}_{n}(t) at ct = 0", \
"fort.2011" u 1:5  w l    ti ""
pause  1
p[0:40][0:0.32] "fort.2501" u 1:2 ti "ct = 2000 fm"
pause 0.4
p[0:40][0:0.32] "fort.2502" u 1:2 ti "ct = 4000 fm"
pause 0.4
p[0:40][0:0.32] "fort.2503" u 1:2 ti "ct = 6000 fm"


