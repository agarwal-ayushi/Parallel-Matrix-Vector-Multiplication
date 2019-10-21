#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <cmath>
#include <omp.h>


#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW(id+1, p, n) - 1 )
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW(id+1, p, n) - BLOCK_LOW(id, p, n))
#define BLOCK_OWNER(index, p, n) (((p) * ((index)+1)-1)/(n))

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
void matrixVectorMulOpenMP(double* result, double* matrix, double* vector, int rowCount, int colCount){
  int i, j;
  #pragma omp parallel shared(matrix, vector, result) private(i,j)
  {
  #pragma omp for schedule(static)
  for (i = 0; i < rowCount; i++) {
		result[i] = 0;
		for (j = 0; j < colCount; j++) {
			//printf("Thread %d did %d,%d\n", omp_get_thread_num(), i, j);
			result[i] += matrix[i*colCount+j]*vector[j];
    }
  }
  }
}
void matrixVectorMulCheckerOpenMP(double* result_openmp, double* matrix, double* vector, int rowCount, int colCount){
  int i, j, id, num_threads, grid[2];
	for (int r = 0; r < rowCount; r++) {
		result_openmp[r] =0.0;
	}
	printMatrix(matrix, rowCount, colCount);
	printVector(vector, colCount);
  #pragma omp parallel shared(matrix, vector, result_openmp, rowCount, colCount) private(i,j, id) firstprivate(grid)
  {
	num_threads = omp_get_num_threads();

	if (rowCount > colCount){
		grid[1] = sqrt(num_threads);
		grid[0] = num_threads/grid[1];
	} else {
		grid[0] = sqrt(num_threads);
		grid[1] = num_threads/grid[0];
	}
	if (omp_get_thread_num() == 0) cout << grid[0] << "\t" << grid[1] << endl;
	#pragma omp barrier
	//int* pNumRows = new int [num_threads];
	//int* pNumCols = new int [num_threads];
	id = omp_get_thread_num();
	int pNumRows = BLOCK_SIZE(id, grid[0], rowCount);
	int pNumCols = BLOCK_SIZE(id, grid[1], colCount);
	int* pRowInd = new int [grid[0]];
	int* pRowNum = new int [grid[0]];
	int* pColInd = new int [grid[1]];
	int* pColNum = new int [grid[1]];
	pRowInd[0] = 0;
	pRowNum[0] = BLOCK_SIZE(0, grid[0], rowCount);
	pColInd[0] = 0;
	pColNum[0] = BLOCK_SIZE(0, grid[1], colCount);
	for (int i = 0; i < grid[0]; i++) {
		pRowInd[i] = pRowInd[i-1] + pRowNum[i-1];
		pRowNum[i] = BLOCK_SIZE(i, grid[0], rowCount);
	}
	for (int i = 0; i < grid[1]; i++) {
		pColInd[i] = pColInd[i-1] + pColNum[i-1];
		pColNum[i] = BLOCK_SIZE(i, grid[1], colCount);
	}
	 if (omp_get_thread_num() == 4) {
		for (int i = 0; i < grid[1]; i++) {
	 		cout << pColInd[i] << "\t" << pColNum[i] << endl;
			//cout << pColInd[id/grid[0]]+i << endl;
	 	}
	 	for (int i = 0; i < grid[0]; i++) {
	 		cout << pRowInd[i] << "\t" << pRowNum[i] << endl ;
			//cout <<pRowInd[id%grid[0]]+i<< endl;
	 	}
	}
	// }
	double* tempResult = new double [pNumRows];
	//printf("%d\n", BLOCK_LOW(id, num_threads/2, colCount));
	//printf("num_threads = %d, thread = %d, rows = %d, cols= %d\n", num_threads, id, pNumRows, pNumCols);
	#pragma omp barrier
	for (i = 0; i < pNumRows; i++) {
		tempResult[i] = 0;
		for (j = 0; j < pNumCols; j++) {
			//printf("Thread %d did row = %d Col= %d\n", id , i, j);
			tempResult[i] += matrix[(pRowInd[id%grid[0]]+i)*colCount + pColInd[id/grid[0]]+j]*vector[pColInd[id/grid[0]]+j];
		}
	}
	if (omp_get_thread_num() == 3) printVector(tempResult, pNumRows);

	#pragma omp barrier
	#pragma omp critical
	for (int i = 0; i < grid[0]; i++){
		if (id%grid[0] == i){
			for (int j = 0; j < pNumRows; j++) {
				result_openmp[pRowInd[i]+j] += tempResult[j];
			}
		}
	}
}
}
void matrixVectorMulCheckerboardOpenMP(double* result, double* matrix, double* vector, int rowCount, int colCount){
	int i, j, num_threads;
	num_threads = omp_get_num_threads();
  #pragma omp parallel shared(matrix, vector, result,num_threads) private(i,j)
  {
  #pragma omp for schedule(static, num_threads) collapse(2)
  for (i = 0; i < rowCount; i++) {
		for (j = 0; j < colCount; j++) {
			printf("Thread %d did %d,%d\n", omp_get_thread_num(), i, j);
			result[i] += matrix[i*colCount+j]*vector[j];
    }
  }
  }
}
void matrixVectorMulRowOpenMP(double* result, double* matrix, double* vector, int rowCount, int colCount){
  int i, j;
  #pragma omp parallel default(none) shared(i, matrix, vector, result, rowCount, colCount) private(j)
  {
  #pragma omp for schedule(static)
  for (i = 0; i < rowCount; i++) {
		printf("Thread %d did row %d\n", omp_get_thread_num(), i);
		result[i] = 0;
		for (j = 0; j < colCount; j++) {
			printf("Thread %d did Col %d\n", omp_get_thread_num(), j);
			result[i] += matrix[i*colCount+j]*vector[j];
    }
  }
  }
}
void matrixVectorMulColOpenMP(double* result, double* matrix, double* vector, int rowCount, int colCount){
  int i, j; double localResult=0.0;
  #pragma omp parallel default(none) shared(j, matrix, vector, result, rowCount, colCount) private(i, localResult)
  {
  for (i = 0; i < rowCount; i++) {
		result[i] = 0;
		printf("Thread %d did row %d\n", omp_get_thread_num(), i);
		#pragma omp for schedule(static)
		for (j = 0; j < colCount; j++) {
			printf("Thread %d did Col %d\n", omp_get_thread_num(), j);
			localResult += matrix[i*colCount+j]*vector[j];
    }
		#pragma omp critical
		{
		result[i] += localResult;
		localResult = 0.0;
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
  //struct timespec start, end;
	double start, end;
	double diff_parallel=0.0, diff_serial=0.0, speedup=0.0;
	//inputInit(rows, cols, size);
	double* matrix = new double [rows*cols];
	double *vector = new double[size];
	double *result = new double[rows];
	double* result_openmp = new double [rows];
  RandomInit(matrix, vector, rows, cols);

	//Serial Matrix Vector multiplication timing

	//clock_gettime(CLOCK_MONOTONIC_RAW, &start);
	start = omp_get_wtime();
	matrixVectorMulSerial(result, matrix, vector, rows, cols);
	//clock_gettime(CLOCK_MONOTONIC_RAW, &end);
	end = omp_get_wtime();
	diff_serial = end - start;
	//printVector(result, rows);
	//diff_serial = double((end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec)/1000000000.0);

	//OpenMP normal multiplication timing
  //clock_gettime(CLOCK_MONOTONIC_RAW, &start);
	start = omp_get_wtime();
	matrixVectorMulCheckerOpenMP(result_openmp, matrix, vector, rows, cols);
	//matrixVectorMulOpenMP(result, matrix, vector, rows, cols);
	//clock_gettime(CLOCK_MONOTONIC_RAW, &end);
	printVector(result_openmp, rows);
	end = omp_get_wtime();
	diff_parallel = end - start;
	//diff_parallel = double((end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec)/1000000000.0);
	testResult(matrix, vector, result_openmp, rows, cols);

	//speedup Serial vs OpenMP
	speedup = diff_serial/diff_parallel;
	printf("rows = %d\tcols = %d\tSerial Time = %fs\tOpenMP Time = %fs\tSpeedup = %f\n", rows, cols, diff_serial, diff_parallel, speedup);
	processTerminate(matrix, vector, result);
  return 0;
}
