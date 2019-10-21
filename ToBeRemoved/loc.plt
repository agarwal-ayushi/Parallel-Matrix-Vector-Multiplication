#!/usr/bin/gnuplot --persist

set title "Sequential vs Parallel (with threads)"
set xlabel "#threads"
set ylabel "Execution Time"
set grid
set style arrow 1 nohead
plot "result_parallel.dat" u (column(0)):2:xtic(1) smooth csplines w l linestyle 1 linewidth 3 title "Parralel Time vs #threads", \
	"" u (column(0)):2:xtic(1) with points linestyle 2 linewidth 3 title "", \
	"result_seq.dat" u (column(0)):2 with lines linewidth 2 lc rgb "red" title "Seq Time"
pause -1


set term x11 0
plot sin(x)
set term x11 1


set term wxt 1
