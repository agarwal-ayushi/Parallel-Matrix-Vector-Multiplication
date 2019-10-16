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
void printMatrix(double** matrix, int rowCount, int colCount) {
	for (int i = 0; i < rowCount; i++){
		for (int j = 0; j < colCount; j++){
			cout << matrix[i][j] << "\t";
		}
		cout << "\n";
	}
	cout << endl;
}
void printVector(double *vector, int size) {
	for (int i = 0; i < size; i++){
		cout << vector[i] << "\t";
	}
	cout << endl;
}
void matrixInit(double** matrix, int rowCount, int colCount){
	for (int i = 0; i < rowCount; i++){
		for (int j = 0; j < colCount; j++){
			matrix[i][j] = i*colCount + j;
		}
	}
}
void RandomInit(double** matrix, double* vector, int rowCount, int colCount){
	unsigned int seed;
	for (int i = 0; i < rowCount; i++){
		for (int j=0; j < colCount; j++){
			matrix[i][j] = double(rand_r(&seed) % 100);
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
void matrixVectorMul(double* result, double** matrix, double* vector, int rowCount, int colCount) {
	for (int i = 0; i < rowCount; i++) {
		result[i] = 0;
		for (int j = 0; j < colCount; j++) {
			result[i] += matrix[i][j]*vector[j];
		}
	}
}
void processTerminate(double** matrix, double* vector, double* result){
	delete [] matrix;
	delete [] vector;
	delete [] result;
}
int main(int argc, char* argv[]) {
	int rows, cols, size;
	rows = stoi(argv[1]);	cols = stoi(argv[2]);	size = stoi(argv[3]);
	if (size != cols) {
		cout << "The input vector is incompatible with the matrix, Please correct input\n";
		exit(1);
	}
	struct timespec start, end;
	double diff=0.0;
	//inputInit(rows, cols, size);
	double **matrix = new double*[rows];
	for (int i = 0; i < rows; i++){
		matrix[i] = new double[cols];
	}
	double *vector = new double[size];
	double *result = new double[rows];

	//matrixInit(matrix, rows, cols);
	//vectorInit(vector, size);
	RandomInit(matrix, vector, rows, cols);
	//printMatrix(matrix, rows, cols);
	//printVector(vector, size);

	clock_gettime(CLOCK_MONOTONIC_RAW, &start);
	matrixVectorMul(result, matrix, vector, rows, cols);
	clock_gettime(CLOCK_MONOTONIC_RAW, &end);
	diff = double((end.tv_sec - start.tv_sec) + (end.tv_nsec - start
	.tv_nsec)/1000000000.0);
	//printVector(result, rows);
	printf("%d\t%d\t%d\t%f\n", rows, cols, size, diff);
	processTerminate(matrix, vector, result);
	return 0;
}
