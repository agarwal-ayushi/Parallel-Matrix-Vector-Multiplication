#!/usr/bin/bash
rm -rf resultsRow_MPI.dat
file="resultsRow_MPI.dat"
exec 3<> $file
#Row-Wise Striped Matrix OpenMP Results
mpicxx mvm_rs_mpi.cpp -o mvm_rs_mpi

echo "\nNP = 1" >&3
mpirun -np 1 ./mvm_rs_mpi 10000 10000 10000 >&3

echo "\nNP =2" >&3
mpirun -np 2 ./mvm_rs_mpi 10000 10000 10000 >&3

echo "\nNP =3" >&3
mpirun -np 3 ./mvm_rs_mpi 10000 10000 10000 >&3

echo "\nNP =4" >&3
mpirun -np 4 ./mvm_rs_mpi 10000 10000 10000 >&3

echo "\nNP =8" >&3
mpirun -np 8 ./mvm_rs_mpi 10000 10000 10000 >&3
