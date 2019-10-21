#!/usr/bin/bash
rm -rf resultsColOpenmp.dat
file="resultsColOpenmp.dat"
exec 3<> $file
#Row-Wise Striped Matrix OpenMP Results
g++ -fopenmp columnStripe_openmp.cpp -o columnStripe_openmp

export OMP_NUM_THREADS=1
echo "\nOMP_NUM_THREADS =1" >&3
./columnStripe_openmp 1000 1000 1000 >&3

echo "\nOMP_NUM_THREADS =2" >&3
export OMP_NUM_THREADS=2
./columnStripe_openmp 1000 1000 1000 >&3

echo "\nOMP_NUM_THREADS =3" >&3
export OMP_NUM_THREADS=3
./columnStripe_openmp 1000 1000 1000 >&3

echo "\nOMP_NUM_THREADS =4" >&3
export OMP_NUM_THREADS=4
./columnStripe_openmp 1000 1000 1000 >&3

echo "\nOMP_NUM_THREADS =8" >&3
export OMP_NUM_THREADS=8
./columnStripe_openmp 1000 1000 1000 >&3
