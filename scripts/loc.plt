#!/usr/bin/gnuplot --persist

set term wxt 0
set title "Speedup with Parallel OpenMP w.r.t Serial"
set xlabel "#omp_threads"
set ylabel "Speedup"
set grid
set style arrow 1 nohead
plot "result_openmp.dat" u (column(0)):3:xtic(1) smooth csplines w l linestyle 1 linewidth 3 title "Row-wise OpenMP", \
	"" u (column(0)):3:xtic(1) with points linestyle 1 linewidth 3 title "", \
	"" u (column(0)):5:xtic(1) smooth csplines w l linestyle 2 linewidth 3 title "Column-wise OpenMP", \
	"" u (column(0)):5:xtic(1) with points linestyle 2 linewidth 3 title "", \
	"" u (column(0)):7:xtic(1) smooth csplines w l linestyle 3 linewidth 3 title "Checkerboard OpenMP", \
	"" u (column(0)):7:xtic(1) with points linestyle 3 linewidth 3 title "",

set term wxt 1
set title "Speedup with Increasing processors of Parallel OpenMP"
set xlabel "#threads"
set ylabel "Execution Time"
set grid
set style arrow 1 nohead
plot "result_openmp.dat" u (column(0)):2:xtic(1) smooth csplines w l linestyle 1 linewidth 3 title "Row-wise OpenMP", \
	"" u (column(0)):2:xtic(1) with points linestyle 1 linewidth 3 title "", \
	"" u (column(0)):4:xtic(1) smooth csplines w l linestyle 2 linewidth 3 title "Column-wise OpenMP", \
	"" u (column(0)):4:xtic(1) with points linestyle 2 linewidth 3 title "", \
	"" u (column(0)):6:xtic(1) smooth csplines w l linestyle 3 linewidth 3 title "Checkerboard OpenMP", \
	"" u (column(0)):6:xtic(1) with points linestyle 3 linewidth 3 title "",

	set term wxt 2
	set title "Speedup with Parallel MPI w.r.t Serial"
	set xlabel "#omp_threads"
	set ylabel "Speedup"
	set grid
	set style arrow 1 nohead
	plot "result_mpi.dat" u (column(0)):3:xtic(1) smooth csplines w l linestyle 1 linewidth 3 title "Row-wise MPI", \
		"" u (column(0)):3:xtic(1) with points linestyle 1 linewidth 3 title "", \
		"" u (column(0)):5:xtic(1) smooth csplines w l linestyle 2 linewidth 3 title "Column-wise MPI", \
		"" u (column(0)):5:xtic(1) with points linestyle 2 linewidth 3 title "", \
		"" u (column(0)):7:xtic(1) smooth csplines w l linestyle 3 linewidth 3 title "Checkerboard MPI", \
		"" u (column(0)):7:xtic(1) with points linestyle 3 linewidth 3 title "",

	set term wxt 3
	set title "Speedup with Increasing processors of Parallel MPI"
	set xlabel "#threads"
	set ylabel "Execution Time"
	set grid
	set style arrow 1 nohead
	plot "result_mpi.dat" u (column(0)):2:xtic(1) smooth csplines w l linestyle 1 linewidth 3 title "Row-wise MPI", \
		"" u (column(0)):2:xtic(1) with points linestyle 1 linewidth 3 title "", \
		"" u (column(0)):4:xtic(1) smooth csplines w l linestyle 2 linewidth 3 title "Column-wise MPI", \
		"" u (column(0)):4:xtic(1) with points linestyle 2 linewidth 3 title "", \
		"" u (column(0)):6:xtic(1) smooth csplines w l linestyle 3 linewidth 3 title "Checkerboard MPI", \
		"" u (column(0)):6:xtic(1) with points linestyle 3 linewidth 3 title "",


pause -1
