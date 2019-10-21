#include <stdio.h>
#include <iostream>
#include <time.h>
#include <stdlib.h>

using namespace std;
void inputInit(int &rows, int &cols, int &size){
	cout << "Enter the size of the matrix and vector in the form of rows, cols\n";
	cout << "Number of rows in Matrix=";
	cin >> rows;
	cout << "Number of columns in the Matrix=";
	cin >> cols;
	cout << "Number of elements in Vector=";
	cin >> size;
	if (size != cols) {
		cout << "The input vector is incompatible with the matrix, Please correct input\n";
		exit(1);
	}
}
void printMatrix(double* matrix, int rowCount, int colCount) {
	for (int i = 0; i < rowCount; i++){
		for (int j = 0; j < colCount; j++){
			printf("%f\t", matrix[i*colCount +j]);
		}
		printf("\n");
	}
	fflush(stdout);
}
void printVector(double *vector, int size) {
	for (int i = 0; i < size; i++){
		printf("%f\t", vector[i]);
	}
	printf("\n");
	fflush(stdout);
}
void matrixInit(double* matrix, int rowCount, int colCount){
	for (int i = 0; i < rowCount; i++){
		for (int j = 0; j < colCount; j++){
			matrix[i*colCount+j] = i*colCount + j;
		}
	}
}
void vectorInit(double* vector, int size){
	for (int i = 0; i < size; i++) {
		vector[i] = i;
	}
}
void RandomInit(double* matrix, double* vector, int rowCount, int colCount){
	unsigned int seed;
	for (int i = 0; i < rowCount; i++){
		for (int j=0; j < colCount; j++){
			matrix[i*colCount+j] = double(rand_r(&seed) % 100);
		}
	}
	for (int i = 0; i < colCount; i++){
		vector[i] = double(rand_r(&seed) % 100);
	}
}
void matrixVectorMulSerial(double* result, double* matrix, double* vector, int rowCount, int colCount) {
	for (int i = 0; i < rowCount; i++) {
		result[i] = 0;
		for (int j = 0; j < colCount; j++) {
			result[i] += matrix[i*colCount+j]*vector[j];
		}
	}
}
void processTerminate(double* matrix, double* vector, double* result){
	delete [] matrix;
	delete [] vector;
	delete [] result;
}
int main(int argc, char* argv[]) {
	int rows, cols, size;
	rows = stoi(argv[1]);	cols = stoi(argv[2]);	size = stoi(argv[3]);
	if (size != cols) {
		printf("The input vector is incompatible with the matrix, Please correct input\n");
		exit(1);
	}
	struct timespec start, end;
	double diff=0.0;
	double *matrix = new double[rows*cols];
	double *vector = new double[size];
	double *result = new double[rows];

	RandomInit(matrix, vector, rows, cols);
	clock_gettime(CLOCK_MONOTONIC_RAW, &start);
	matrixVectorMulSerial(result, matrix, vector, rows, cols);
	clock_gettime(CLOCK_MONOTONIC_RAW, &end);
	diff = double((end.tv_sec - start.tv_sec) + (end.tv_nsec - start
	.tv_nsec)/1000000000.0);
	printf("rows = %d\tcols = %d\tTime = %fs\n", rows, cols, diff);
	processTerminate(matrix, vector, result);
	return 0;
}
