#!/bin/bash

rm -f result_openmp.dat result_mpi.dat
sed -n -e '/Row/,/Column/ p' results_openmp.txt | grep -E "OMP_NUM_THREADS" > temp.txt
touch "loc1.dat"
touch "loc2.dat"
touch "loc3.dat"
touch "loc4.dat"
touch "loc5.dat"
touch "loc6.dat"
touch "loc7.dat"
echo "#THREADS	" "Row Time		" "Row Speedup		" "Col Time		" "Col Speedup		" "Checker Time	" "Checker Speedup	" > result_openmp.dat

set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc1.dat
	i=$((i+1))
done < temp.txt

sed -n -e '/Row/,/Column/ p' results_openmp.txt | grep "OpenMP Time" > temp_row.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc2.dat
	i=$((i+1))
done < temp_row.txt

sed -n -e '/Row/,/Column/ p' results_openmp.txt | grep "Speedup" > temp_row.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc5.dat
	i=$((i+1))
done < temp_row.txt

sed -n -e '/Column/,/Checkerboard/ p' results_openmp.txt | grep "OpenMP Time" > temp_col.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc3.dat
	i=$((i+1))
done < temp_col.txt

sed -n -e '/Column/,/Checkerboard/ p' results_openmp.txt | grep "Speedup" > temp_col.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc6.dat
	i=$((i+1))
done < temp_col.txt

sed -n -e '/Checkerboard/,// p' results_openmp.txt | grep "OpenMP Time" > temp_check.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc4.dat
	i=$((i+1))
done < temp_check.txt

sed -n -e '/Checkerboard/,// p' results_openmp.txt | grep "Speedup" > temp_check.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc7.dat
	i=$((i+1))
done < temp_check.txt

paste <(awk '{print $1"	"}' loc1.dat ) <(awk '{print $1"	"}' loc2.dat ) <(awk '{print $1"	" }' loc5.dat ) <(awk '{print $1"	" }' loc3.dat ) <(awk '{print $1"	" }' loc6.dat ) <(awk '{print $1}' loc4.dat ) <(awk '{print $1"		" }' loc7.dat ) >>result_openmp.dat
rm -f loc5.dat loc6.dat loc7.dat loc3.dat loc1.dat loc2.dat loc4.dat temp.txt temp_row.txt temp_col.txt temp_check.txt

# Print MPI result
sed -n -e '/Row/,/Column/ p' results_mpi.txt | grep -E "NP" > temp.txt
touch "loc1.dat"
touch "loc2.dat"
touch "loc3.dat"
touch "loc4.dat"
touch "loc5.dat"
touch "loc6.dat"
touch "loc7.dat"
echo "#NP	" "Row Time		" "Row Speedup		" "Col Time		" "Col Speedup		" "Checker Time	" "Checker Speedup	" > result_mpi.dat

set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc1.dat
	i=$((i+1))
done < temp.txt

sed -n -e '/Row/,/Column/ p' results_mpi.txt | grep "MPI Time" > temp_row.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc2.dat
	i=$((i+1))
done < temp_row.txt

sed -n -e '/Row/,/Column/ p' results_mpi.txt | grep "Speedup" > temp_row.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc5.dat
	i=$((i+1))
done < temp_row.txt

sed -n -e '/Column/,/Checkerboard/ p' results_mpi.txt | grep "MPI Time" > temp_col.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc3.dat
	i=$((i+1))
done < temp_col.txt

sed -n -e '/Column/,/Checkerboard/ p' results_mpi.txt | grep "Speedup" > temp_col.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc6.dat
	i=$((i+1))
done < temp_col.txt

sed -n -e '/Checkerboard/,/\\n/ p' results_mpi.txt | grep "MPI Time" > temp_check.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc4.dat
	i=$((i+1))
done < temp_check.txt

sed -n -e '/Checkerboard/,/\\n/ p' results_mpi.txt | grep "Speedup" > temp_check.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc7.dat
	i=$((i+1))
done < temp_check.txt

paste <(awk '{print $1"	"}' loc1.dat ) <(awk '{print $1"	"}' loc2.dat ) <(awk '{print $1"	" }' loc5.dat ) <(awk '{print $1"	" }' loc3.dat ) <(awk '{print $1"	" }' loc6.dat ) <(awk '{print $1}' loc4.dat ) <(awk '{print $1"		" }' loc7.dat ) >>result_mpi.dat
rm -f loc5.dat loc6.dat loc7.dat loc3.dat loc1.dat loc2.dat loc4.dat temp.txt temp_row.txt temp_col.txt temp_check.txt
gnuplot -c scripts/loc.plt
