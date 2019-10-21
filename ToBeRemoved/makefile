CC = gcc
CXX = g++

CFLAGS = -g -Wall

all: pi pi_mt report1_new.txt

pi: pi.cpp
	$(CXX) $(CFLAGS) -o pi pi.cpp
pi_mt: pi_mt.cpp
	$(CXX) $(CFLAGS) pi_mt.cpp -o pi_mt -lpthread 

clean: 
	rm pi pi_mt

