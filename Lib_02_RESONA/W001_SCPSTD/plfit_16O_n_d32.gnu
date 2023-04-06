
set rmargin 0.3
set tmargin 0.3
set size 0.8, 0.8

set xlabel 'E  (MeV)'
set ylabel 'd{/Symbol s}(E)/dE'

set style line 1 lt 1 lw 2 lc rgb 'black'
set style line 2 lt 2 lw 2 lc rgb 'red'
set style line 3 lt 3 lw 2 lc rgb 'blue'
set style line 4 lt 4 lw 2 lc rgb 'purple'
set style line 5 lt 5 lw 2 lc rgb 'magenta'
set style line 6 lt 6 lw 2 lc rgb 'cyan'

e1 = 0.415563
g1 = 0.0166501
###c1 = 0.01
###d1 = 0.01
ea = 0.36
eb = 0.48

f1(x)= (g1/2) / ( (g1/2)**2 + (x-e1)**2 ) ###+ c1*x + d1
fit[ea : eb]  f1(x) 'fort.290023'  u   1:3 via e1, g1

p[ea : eb][0:144] \
'fort.290023' u 1:3 ti '^{16}O + neutron in d_{3/2}', \
f1(x) w l ls 3 ti 'fitting'

pause 2

###set title "Scattering of ^{16}O + p"
set terminal postscript eps enhanced color font "Century, 20"
#set output "hogee.eps"
set output "| epstopdf -f -o=hogee.pdf"
replot
unset output
set terminal x11
unset label

reset

