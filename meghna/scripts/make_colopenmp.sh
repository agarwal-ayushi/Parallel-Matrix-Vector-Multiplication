#!/usr/bin/bash
rm -rf resultsColOpenmp_lab.dat
file="resultsColOpenmp_lab.dat"
exec 3<> $file
#Col-Wise Striped Matrix OpenMP Results

export OMP_NUM_THREADS=1
echo "\nOMP_NUM_THREADS =1" >&3
./columnStripe_openmp 10000 10000 10000 >&3

echo "\nOMP_NUM_THREADS =2" >&3
export OMP_NUM_THREADS=2
./columnStripe_openmp 10000 10000 10000 >&3

echo "\nOMP_NUM_THREADS =3" >&3
export OMP_NUM_THREADS=3
./columnStripe_openmp 10000 10000 10000 >&3

echo "\nOMP_NUM_THREADS =4" >&3
export OMP_NUM_THREADS=4
./columnStripe_openmp 10000 10000 10000 >&3

echo "\nOMP_NUM_THREADS =8" >&3
export OMP_NUM_THREADS=8
./columnStripe_openmp 10000 10000 10000 >&3

echo "\nOMP_NUM_THREADS =16" >&3
export OMP_NUM_THREADS=16
./columnStripe_openmp 10000 10000 10000 >&3
