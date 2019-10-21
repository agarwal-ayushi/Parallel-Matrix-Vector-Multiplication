# Open-MPI-Matrix-Vector-Multiplication

This repository contains the parallel Open MPI and OpenMP implementation of Matrix Vector Multiplication using three methods:
1. Row-wise striped
2. Column-Wise Striped
3. Checkerboard Striped

To run, please do the following:

`sh make.sh`
- This will first make all the files present in the `src/`
- Then a report is generated for OpenMP and MPI separately
- Then you can plot graph for OpenMP and MPI speedup obtained vs Serial implementation and as the number of processes is increased
