#!/bin/bash
size=$cols
if [ $FUNCTIONAL_MODE -eq 0 ];
 then
  sh scripts/report.sh
  ./scripts/plot.sh
elif [ $FUNCTIONAL_MODE -eq 1 ]
  then
  echo "Performing Row-Wise Matrix-Vector Multiplication with OpenMP and MPI\n"
  export OMP_NUM_THREADS=$proc
  ./rowStripe_openmp $rows $cols $size
  mpirun -np $proc ./mvm_rs_mpi $rows $cols $size
elif [ $FUNCTIONAL_MODE -eq 2 ]
  then
  echo "Performing Column-Wise Matrix-Vector Multiplication with OpenMP and MPI\n"
  export OMP_NUM_THREADS=$proc
  ./columnStripe_openmp $rows $cols $size
  mpirun -np $proc ./mvm_cs_mpi $rows $cols $size
elif [ $FUNCTIONAL_MODE -eq 3 ]
  then
  echo "Performing Chckerboard Matrix-Vector Multiplication with OpenMP and MPI\n"
  export OMP_NUM_THREADS=$proc
  ./checkerboard_openmp $rows $cols $size $OMP_NUM_THREADS
  mpirun -np $proc ./mvm_checkerboard $rows $cols $size
fi
