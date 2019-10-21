#!/usr/bin/bash
rm -rf resultsRowOpenmp.dat
file="resultsRowOpenmp.dat"
exec 3<> $file
#Row-Wise Striped Matrix OpenMP Results
g++ -fopenmp rowStripe_openmp.cpp -o rowStripe_openmp

export OMP_NUM_THREADS=1
echo "\nOMP_NUM_THREADS =1" >&3
./rowStripe_openmp 1000 1000 1000 >&3
echo "\nOMP_NUM_THREADS =2" >&3
export OMP_NUM_THREADS=2
./rowStripe_openmp 1000 1000 1000 >&3
echo "\nOMP_NUM_THREADS =3" >&3
export OMP_NUM_THREADS=3
./rowStripe_openmp 1000 1000 1000 >&3
echo "\nOMP_NUM_THREADS =4" >&3
export OMP_NUM_THREADS=4
./rowStripe_openmp 1000 1000 1000 >&3
echo "\nOMP_NUM_THREADS =8" >&3
export OMP_NUM_THREADS=8
./rowStripe_openmp 1000 1000 1000 >&3
