#!/usr/bin/bash
rm -rf resultsCol_MPI.dat
file="resultsCol_MPI.dat"
exec 3<> $file
#Column-Wise Striped Matrix OpenMP Results
mpicxx mvm_cs_mpi.cpp -o mvm_cs_mpi

echo "\nNP = 1" >&3
mpirun -np 1 ./mvm_cs_mpi 10000 10000 10000 >&3

echo "\nNP =2" >&3
mpirun -np 2 ./mvm_cs_mpi 10000 10000 10000 >&3

echo "\nNP =3" >&3
mpirun -np 3 ./mvm_cs_mpi 10000 10000 10000 >&3

echo "\nNP =4" >&3
mpirun -np 4 ./mvm_cs_mpi 10000 10000 10000 >&3

echo "\nNP =8" >&3
mpirun -np 8 ./mvm_cs_mpi 10000 10000 10000 >&3
