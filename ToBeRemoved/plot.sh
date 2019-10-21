#!/bin/bash

rm -f loc3.dat result_seq.dat result_parallel.dat loc1.dat loc2.dat temp.txt temp1.txt
grep -E "threads" report1_new.txt > temp.txt
touch "loc1.dat"
touch "loc2.dat"
touch "loc3.dat"
echo "#THREADS     " "Exec Time" > result_parallel.dat
echo "#THREADS     " "Exec Time" > result_seq.dat

set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc1.dat
	i=$((i+1))
done < temp.txt
 
grep -E "time" report1_new.txt > temp1.txt
set i = 1
IFS="="
while read line; do
	read -ra ADDR <<< "$line"
	echo ${ADDR[1]} >> loc2.dat
	i=$((i+1))
done < temp1.txt

for i in {1..9}; do sed -n '1p' loc2.dat >> loc3.dat; done

sed -i".bak" '1d' loc2.dat 
paste <(awk '{print $1}' loc1.dat ) <(awk '{print $1}' loc2.dat ) >>result_parallel.dat
paste <(awk '{print $1}' loc1.dat ) <(awk '{print $1}' loc3.dat ) >>result_seq.dat

rm -f temp.txt temp1.txt loc1.dat loc2.dat loc3.dat loc2.dat.bak
gnuplot -c loc.plt

