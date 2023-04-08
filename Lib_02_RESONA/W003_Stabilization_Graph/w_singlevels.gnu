set tmargin 0.3
set rmargin 0.3


p[][0:]  \
"Fine.dat"  u  1:2 w  p lt 1 ti "^{16}O + neut.(d_{3/2}),", \
"Fine.dat"  u  1:3 w  p lt 1 ti "", \
"Fine.dat"  u  1:4 w  p lt 6 lc "black"   ti "", \
"Fine.dat"  u  1:5 w  p lt 1 ti "", \
"Fine.dat"  u  1:6 w  p lt 1 ti ""
pause 1

set size 0.6, 0.6
set xlabel "R_{max}  [fm]"
set ylabel "E_{SP}  [MeV]"
set terminal postscript eps enhanced color font "Century, 20"
set output "| epstopdf -f -o=hogee.pdf"
replot
unset output
set terminal x11
unset label
