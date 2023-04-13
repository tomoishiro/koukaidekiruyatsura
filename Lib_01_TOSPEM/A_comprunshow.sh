#
rm  a.out
gfortran  TOSPEM.f90
./a.out
cat RESULT_RUN*.dat
#
#See "PARAM.inp" for input parameters.
#
