#!/bin/bash
make clean
make all
rm -f report1_new.txt
{
echo "PI converges in N = 1000000\n"

echo "Output of Sequential Code\n"
./pi 1000000

echo "\n"

echo "Output of Parallel Code\n"
./pi_mt 1000000 1
echo "\n"
./pi_mt 1000000 2
echo "\n"
./pi_mt 1000000 4
echo "\n"
./pi_mt 1000000 6
echo "\n"
./pi_mt 1000000 8
echo "\n"
./pi_mt 1000000 16
echo "\n"
./pi_mt 1000000 32
echo "\n"
./pi_mt 1000000 40
echo "\n"
./pi_mt 1000000 100

} > report1_new.txt
