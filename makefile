CC = gcc
CXX = g++
MPICXX = mpicxx
CFLAGS = -g -Wall
CFLAGS_opemp = -g -Wall -fopenmp

all: serial rowStripe_openmp columnStripe_openmp checkerboard_openmp mvm_rs_mpi mvm_cs_mpi mvm_checkerboard

serial: serial.cpp
	$(CXX) $(CFLAGS) serial.cpp -o serial
rowStripe_openmp: rowStripe_openmp.cpp
	$(CXX) $(CFLAGS_opemp) rowStripe_openmp.cpp -o rowStripe_openmp
columnStripe_openmp: columnStripe_openmp.cpp
	$(CXX) $(CFLAGS_opemp) columnStripe_openmp.cpp -o columnStripe_openmp
checkerboard_openmp: checkerboard_openmp.cpp
	$(CXX) $(CFLAGS_opemp) checkerboard_openmp.cpp -o checkerboard_openmp
mvm_rs_mpi: mvm_rs_mpi.cpp
	$(MPICXX) mvm_rs_mpi.cpp -o mvm_rs_mpi
mvm_cs_mpi: mvm_cs_mpi.cpp
	$(MPICXX) mvm_cs_mpi.cpp -o mvm_cs_mpi
mvm_checkerboard: mvm_checkerboard.cpp
	$(MPICXX) mvm_checkerboard.cpp -o mvm_checkerboard
clean:
	rm rowStripe_openmp columnStripe_openmp checkerboard_openmp mvm_rs_mpi mvm_cs_mpi mvm_checkerboard
