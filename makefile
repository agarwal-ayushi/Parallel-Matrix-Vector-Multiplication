CC = gcc
CXX = g++
MPICXX = mpicxx
CFLAGS = -g -Wall
CFLAGS_opemp = -g -Wall -fopenmp

all: serial rowStripe_openmp columnStripe_openmp checkerboard_openmp mvm_rs_mpi mvm_cs_mpi mvm_checkerboard

serial: src/serial.cpp
	$(CXX) $(CFLAGS) src/serial.cpp -o serial
rowStripe_openmp: src/rowStripe_openmp.cpp
	$(CXX) $(CFLAGS_opemp) src/rowStripe_openmp.cpp -o rowStripe_openmp
columnStripe_openmp: src/columnStripe_openmp.cpp
	$(CXX) $(CFLAGS_opemp) src/columnStripe_openmp.cpp -o columnStripe_openmp
checkerboard_openmp: src/checkerboard_openmp.cpp
	$(CXX) $(CFLAGS_opemp) src/checkerboard_openmp.cpp -o checkerboard_openmp
mvm_rs_mpi: src/mvm_rs_mpi.cpp
	$(MPICXX) src/mvm_rs_mpi.cpp -o mvm_rs_mpi
mvm_cs_mpi: src/mvm_cs_mpi.cpp
	$(MPICXX) src/mvm_cs_mpi.cpp -o mvm_cs_mpi
mvm_checkerboard: src/mvm_checkerboard.cpp
	$(MPICXX) src/mvm_checkerboard.cpp -o mvm_checkerboard
clean:
	rm serial rowStripe_openmp columnStripe_openmp checkerboard_openmp mvm_rs_mpi mvm_cs_mpi mvm_checkerboard
