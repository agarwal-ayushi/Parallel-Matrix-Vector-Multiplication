#!/usr/bin/bash
rm -rf resultsChecker_MPI_lab.dat
file="resultsChecker_MPI_lab.dat"
exec 3<> $file
#Column-Wise Striped Matrix OpenMP Results

echo "\nNP =2" >&3
mpirun -np 2 ./mvm_checkerboard 10000 10000 10000 >&3

echo "\nNP =3" >&3
mpirun -np 3 ./mvm_checkerboard 10000 10000 10000 >&3

echo "\nNP =4" >&3
mpirun -np 4 ./mvm_checkerboard 10000 10000 10000 >&3

echo "\nNP =8" >&3
mpirun -np 8 ./mvm_checkerboard 10000 10000 10000 >&3

echo "\nNP =16" >&3
mpirun -np 16 ./mvm_checkerboard 10000 10000 10000 >&3
