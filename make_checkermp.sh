#!/usr/bin/bash
rm -rf resultsCheckerOpenmp.dat
file="resultsCheckerOpenmp.dat"
exec 3<> $file
#Row-Wise Striped Matrix OpenMP Results
g++ -fopenmp checkerboard_openmp.cpp -o checkerboard_openmp

export OMP_NUM_THREADS=1
echo "\nOMP_NUM_THREADS =1" >&3
./checkerboard_openmp 10000 10000 10000 $OMP_NUM_THREADS >&3

echo "\nOMP_NUM_THREADS =2" >&3
export OMP_NUM_THREADS=2
./checkerboard_openmp 10000 10000 10000 $OMP_NUM_THREADS >&3

echo "\nOMP_NUM_THREADS =3" >&3
export OMP_NUM_THREADS=3
./checkerboard_openmp 10000 10000 10000 $OMP_NUM_THREADS >&3

echo "\nOMP_NUM_THREADS =4" >&3
export OMP_NUM_THREADS=4
./checkerboard_openmp 10000 10000 10000 $OMP_NUM_THREADS >&3

echo "\nOMP_NUM_THREADS =8" >&3
export OMP_NUM_THREADS=8
./checkerboard_openmp 10000 10000 10000 $OMP_NUM_THREADS >&3
