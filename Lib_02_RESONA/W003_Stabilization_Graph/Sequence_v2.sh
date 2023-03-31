#!/bin/sh
l=3000000
n=560
z=20
rm  OUTRU*.txt  Rmax_and_En*.txt
while [ $n -ne 1000 ]
do
  echo "(N=${n})-operation"
  sed "4 s/=.../=$n/" SCPSTD.f90 > Main.f90
  gfortran  Main.f90  -o  run.out
  b=`expr $n + $l`
  time ./run.out > OUTRUN_${b}.txt
  mv  Rmax_and_Energies.txt  Rmax_and_Enes_${b}.txt
  n=`expr $n + $z`
  echo "----------------------------------------------------------"
done

#grep  "fm, MeV"  ./Rmax*.txt  >  Dummy.dat
#sed  "s/..Rmax_and_Enes_.....txt:/ /"  Dummy.dat  >  Fine.dat
cat  Rmax_and_Enes_*.txt  >  Fine.dat

