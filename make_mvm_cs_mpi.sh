#!/usr/bin/bash
rm -rf resultsCol_MPI_lab.dat
file="resultsCol_MPI_lab.dat"
exec 3<> $file
#Column-Wise Striped Matrix OpenMP Results

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

echo "\nNP =16" >&3
mpirun -np 16 ./mvm_cs_mpi 10000 10000 10000 >&3
