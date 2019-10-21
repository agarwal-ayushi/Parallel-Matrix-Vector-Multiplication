#!/usr/bin/bash
rm -rf resultsChecker_MPI.dat
file="resultsChecker_MPI.dat"
exec 3<> $file
#Column-Wise Striped Matrix OpenMP Results
mpicxx mvm_checkerboard.cpp -o mvm_checkerboard

echo "\nNP =2" >&3
mpirun -np 2 ./mvm_checkerboard 1000 1000 1000 >&3

echo "\nNP =3" >&3
mpirun -np 3 ./mvm_checkerboard 1000 1000 1000 >&3

echo "\nNP =4" >&3
mpirun -np 4 ./mvm_checkerboard 1000 1000 1000 >&3

echo "\nNP =8" >&3
mpirun -np 8 ./mvm_checkerboard 1000 1000 1000 >&3
