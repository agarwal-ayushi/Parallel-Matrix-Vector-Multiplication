#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <cmath>
#include <omp.h>

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
			printf("%f\t", matrix[i*colCount+j]);
		}
		printf("\n");
	}
	fflush(stdout);
}
void printVector(double* vector, int size) {
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
void vectorInit(double* vector, int size){
	for (int i = 0; i < size; i++) {
		vector[i] = i;
	}
}
void processTerminate(double* matrix, double* vector, double* result){
	delete [] matrix;
	delete [] vector;
	delete [] result;
}
void matrixVectorMulRowOpenMP(double* result, double* matrix, double* vector, int rowCount, int colCount){
  int i, j;
	for (i = 0; i < rowCount; i++) {
		result[i] = 0;
	}
  #pragma omp parallel default(none) shared(i, matrix, vector, result, rowCount, colCount) private(j)
  {
  #pragma omp for schedule(static)
  for (i = 0; i < rowCount; i++) {
		for (j = 0; j < colCount; j++) {
			result[i] += matrix[i*colCount+j]*vector[j];
    }
  }
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
void testResult(double* matrix, double* vector, double* result, int rows, int cols) {
	int equal = 0;
	double* serialResult = new double [rows];
	matrixVectorMulSerial(serialResult, matrix, vector, rows, cols);
	for (int i = 0; i < rows; i++) {
		if (result[i] != serialResult[i]) equal = 1;
	}
	if (equal == 1) printf("Results of parallel and Serial algorithm are not identical. Please check code\n");
	else printf("The result of Serial and Parallel Algorithm are the same. Good Job! \n ");
	delete [] serialResult;
}
int main(int argc, char* argv[]){
  int rows, cols, size;
	rows = stoi(argv[1]);	cols = stoi(argv[2]);	size = stoi(argv[3]);
	if (size != cols) {
		cout << "The input vector is incompatible with the matrix, Please correct input\n";
		exit(1);
	}
	struct timespec start_serial, end_serial;
	double start, end;
	double diff_parallel=0.0, diff_serial=0.0, speedup=0.0, average = 0.0;
	double* matrix = new double [rows*cols];
	double *vector = new double[size];
	double *result = new double[rows];
  RandomInit(matrix, vector, rows, cols);

	//Serial Matrix Vector multiplication timing

	clock_gettime(CLOCK_MONOTONIC_RAW, &start_serial);
	matrixVectorMulSerial(result, matrix, vector, rows, cols);
	clock_gettime(CLOCK_MONOTONIC_RAW, &end_serial);

	diff_serial = double((end_serial.tv_sec - start_serial.tv_sec) + (end_serial.tv_nsec - start_serial.tv_nsec)/1000000000.0);
	//OpenMP normal multiplication timing
	for (int i = 0; i < 100; i++) {
		diff_parallel = 0.0;
		start = omp_get_wtime();
		matrixVectorMulRowOpenMP(result,matrix, vector, rows, cols);
		end = omp_get_wtime();
		diff_parallel = end - start;
		average += diff_parallel;
	}
	average /= 100;
	testResult(matrix, vector, result, rows, cols);

	//speedup Serial vs OpenMP
	speedup = diff_serial/average;
	printf("rows = %d\ncols = %d\nSerial Time = %f\nOpenMP Time = %f\nSpeedup = %f\n", rows, cols, diff_serial, average, speedup);
	//printf("%f\n", diff_parallel);
	processTerminate(matrix, vector, result);
  return 0;
}
