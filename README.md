# Open-MPI-Matrix-Vector-Multiplication

This repository contains the parallel Open MPI and OpenMP implementation of Matrix Vector Multiplication using three methods:
1. Row-wise striped
2. Column-Wise Striped
3. Checkerboard Striped

To run, please do the following:

Please set the following ENV variables on the terminal where you would be running the script.
If you're using bash shell:
1. export rows=\<number of matrix rows\>
2. export cols=\<number of matrix cols\>
3. export proc=\<number of processes\>
4. export FUNCTIONAL_MODE=\<As Below\>

FUNCTIONAL_MODE

	0 => Output the graphs for different number of processes (which can be added in `scripts/report.sh` file)
	
	1 => Use Row-wise method to calculate matrix-vector Multiplication
	
	2 => Use Column-wise method to calculate matrix-vector Multiplication
	
	3 => Use Checkerboard Method to calculate matrix-vector Multiplication

`sh make.sh`
- This will first make all the files present in the `src/`
- The results would be either be in the form of graphs or a report according to the FUNCTIONAL_MODE seleted. 

### Each code also contains a validation step where the output of both parallel and serial implementation of MVM are compared. An error is reported if the results don't match. 
