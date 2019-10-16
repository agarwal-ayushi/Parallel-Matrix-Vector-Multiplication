#!/usr/bin/sh
rm -rf results.dat
file="results.dat"
exec 3<> $file
g++ serial.cpp -o serial
g++ -fopenmp rowStripe_openmp.cpp -o rowStripe_openmp

./serial 100 100 100 >&3
