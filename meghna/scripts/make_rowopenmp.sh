#!/usr/bin/bash
rm -rf resultsRowOpenmp_lab.dat
file="resultsRowOpenmp_lab.dat"
exec 3<> $file
#Row-Wise Striped Matrix OpenMP Results

export OMP_NUM_THREADS=1
echo "\nOMP_NUM_THREADS =1" >&3
./rowStripe_openmp 10000 10000 10000 >&3
echo "\nOMP_NUM_THREADS =2" >&3
export OMP_NUM_THREADS=2
./rowStripe_openmp 10000 10000 10000 >&3
echo "\nOMP_NUM_THREADS =3" >&3
export OMP_NUM_THREADS=3
./rowStripe_openmp 10000 10000 10000 >&3
echo "\nOMP_NUM_THREADS =4" >&3
export OMP_NUM_THREADS=4
./rowStripe_openmp 10000 10000 10000 >&3
echo "\nOMP_NUM_THREADS =8" >&3
export OMP_NUM_THREADS=8
./rowStripe_openmp 10000 10000 10000 >&3
echo "\nOMP_NUM_THREADS =16" >&3
export OMP_NUM_THREADS=16
./rowStripe_openmp 10000 10000 10000 >&3
