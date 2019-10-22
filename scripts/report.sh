#!/bin/bash
make clean
make all
# rows=10000
# cols=10000
# size=10000
size=$cols
# $rows and $cols and $proc would be from ENV variables.
# Please make sure you've set the ENV variables on the terminal running the scripts
rm -rf results_openmp.txt
rm -rf results_mpi.txt
{
  echo "Matrix Size = $rows x $cols "
  echo "Vector Size = $size x 1"
  echo "Row-Wise Striped Matrix-Vector Multiplication with OpenMP\n"
  export OMP_NUM_THREADS=1
  echo "\nOMP_NUM_THREADS =1"
  ./rowStripe_openmp $rows $cols $size
  echo "\nOMP_NUM_THREADS =2"
  export OMP_NUM_THREADS=2
  ./rowStripe_openmp $rows $cols $size
  echo "\nOMP_NUM_THREADS =3"
  export OMP_NUM_THREADS=3
  ./rowStripe_openmp $rows $cols $size
  echo "\nOMP_NUM_THREADS =4"
  export OMP_NUM_THREADS=4
  ./rowStripe_openmp $rows $cols $size
  echo "\nOMP_NUM_THREADS =8"
  export OMP_NUM_THREADS=8
  ./rowStripe_openmp $rows $cols $size
  echo "\nOMP_NUM_THREADS =16"
  export OMP_NUM_THREADS=16
  ./rowStripe_openmp $rows $cols $size

  echo "\nColumn-Wise Striped Matrix-Vector Multiplication with OpenMP\n"
  export OMP_NUM_THREADS=1
  echo "\nOMP_NUM_THREADS =1"
  ./columnStripe_openmp $rows $cols $size
  echo "\nOMP_NUM_THREADS =2"
  export OMP_NUM_THREADS=2
  ./columnStripe_openmp $rows $cols $size
  echo "\nOMP_NUM_THREADS =3"
  export OMP_NUM_THREADS=3
  ./columnStripe_openmp $rows $cols $size
  echo "\nOMP_NUM_THREADS =4"
  export OMP_NUM_THREADS=4
  ./columnStripe_openmp $rows $cols $size
  echo "\nOMP_NUM_THREADS =8"
  export OMP_NUM_THREADS=8
  ./columnStripe_openmp $rows $cols $size
  echo "\nOMP_NUM_THREADS =16"
  export OMP_NUM_THREADS=16
  ./columnStripe_openmp $rows $cols $size

  echo "\nCheckerboard Matrix-Vector Multiplication with OpenMP\n "
  export OMP_NUM_THREADS=1
  echo "\nOMP_NUM_THREADS =1"
  ./checkerboard_openmp $rows $cols $size $OMP_NUM_THREADS
  echo "\nOMP_NUM_THREADS =2"
  export OMP_NUM_THREADS=2
  ./checkerboard_openmp $rows $cols $size $OMP_NUM_THREADS
  echo "\nOMP_NUM_THREADS =3"
  export OMP_NUM_THREADS=3
  ./checkerboard_openmp $rows $cols $size $OMP_NUM_THREADS
  echo "\nOMP_NUM_THREADS =4"
  export OMP_NUM_THREADS=4
  ./checkerboard_openmp $rows $cols $size $OMP_NUM_THREADS
  echo "\nOMP_NUM_THREADS =8"
  export OMP_NUM_THREADS=8
  ./checkerboard_openmp $rows $cols $size $OMP_NUM_THREADS
  echo "\nOMP_NUM_THREADS =16"
  export OMP_NUM_THREADS=16
  ./checkerboard_openmp $rows $cols $size $OMP_NUM_THREADS

} > results_openmp.txt
{
  echo "Matrix Size = $rows x $cols "
  echo "Vector Size = $size x 1"

  echo "Row-Wise Striped Matrix-Vector Multiplication with MPI\n"
  echo "\nNP = 1"
  mpirun -np 1 ./mvm_rs_mpi $rows $cols $size

  echo "\nNP =2"
  mpirun -np 2 ./mvm_rs_mpi $rows $cols $size

  echo "\nNP =3"
  mpirun -np 3 ./mvm_rs_mpi $rows $cols $size

  echo "\nNP =4"
  mpirun -np 4 ./mvm_rs_mpi $rows $cols $size

  echo "\nNP =8"
  mpirun -np 8 ./mvm_rs_mpi $rows $cols $size

  echo "\nNP =16"
  mpirun -np 16 ./mvm_rs_mpi $rows $cols $size

  echo "\nColumn-Wise Striped Matrix-Vector Multiplication with MPI\n"
  echo "\nNP = 1"
  mpirun -np 1 ./mvm_cs_mpi $rows $cols $size

  echo "\nNP =2"
  mpirun -np 2 ./mvm_cs_mpi $rows $cols $size

  echo "\nNP =3"
  mpirun -np 3 ./mvm_cs_mpi $rows $cols $size

  echo "\nNP =4"
  mpirun -np 4 ./mvm_cs_mpi $rows $cols $size

  echo "\nNP =8"
  mpirun -np 8 ./mvm_cs_mpi $rows $cols $size

  echo "\nNP =16"
  mpirun -np 16 ./mvm_cs_mpi $rows $cols $size

  echo "\nCheckerboard Matrix-Vector Multiplication with MPI\n"
  echo "\nNP =1"
  mpirun -np 1 ./mvm_checkerboard $rows $cols $size

  echo "\nNP =2"
  mpirun -np 2 ./mvm_checkerboard $rows $cols $size

  echo "\nNP =3"
  mpirun -np 3 ./mvm_checkerboard $rows $cols $size

  echo "\nNP =4"
  mpirun -np 4 ./mvm_checkerboard $rows $cols $size

  echo "\nNP =8"
  mpirun -np 8 ./mvm_checkerboard $rows $cols $size

  echo "\nNP =16"
  mpirun -np 16 ./mvm_checkerboard $rows $cols $size

} > results_mpi.txt
