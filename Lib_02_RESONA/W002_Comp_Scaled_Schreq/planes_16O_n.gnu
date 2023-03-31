
p[0:12][:3] \
"W100_Poten.dat"          u 1:2, \
"./For_SCPSTD/fort.2001"  u 1:5 w l


pause  1
splot  "W621_ERROR.dat"


pause  1
p[][] \
"W201_WF.dat"     u 1:3 w l, \
"W201_WF.dat"     u 1:4 w l, \
"W301_WF.dat"     u 1:3 w p, \
"W301_WF.dat"     u 1:4 w p

pause  1
p[][] \
"W202_WFDERI.dat" u 1:3 w l, \
"W202_WFDERI.dat" u 1:4 w l, \
"W302_WFDERI.dat" u 1:3 w p, \
"W302_WFDERI.dat" u 1:4 w p

#########################################################################
#########################################################################
#########################################################################
#########################################################################


set terminal postscript eps enhanced color font "Century, 18"
#set terminal postscript eps enhanced monochrome font "Century, 18"
#set output "hogee.eps"
set output "| epstopdf -f -o=hogee3.pdf"


set xlabel "Re(E)  [MeV]"
set xtics -50,0.0001,50
set ylabel "Im(E)  [MeV]"
set ytics -50,0.0001,50

set size 0.7, 0.7
set size ratio 1.0
set pm3d map
set format cb "%2.1t x 10^{%L}"
set isosamples 40 #Use the same pixels for x&y.

dc=0.1
set cbrange[0:dc*10]
set palette defined (0 "white", \
dc*1 "blue", \
dc*2 "turquoise", \
dc*3 "green", \
dc*4 "yellow", \
dc*5 "orange", \
dc*6 "salmon", \
dc*7 "red", \
dc*8 "violet", \
dc*9 "magenta", \
dc*10 "purple")

set title "Matching error (x 6000)"
splot "W621_ERROR.dat" u 1:2:(6000*$5) ti  ""

pause  1

unset output
set terminal x11
unset label
reset

