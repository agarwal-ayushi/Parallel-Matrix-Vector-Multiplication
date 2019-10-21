#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <time.h>
#include <omp.h>
#include <mpi.h>

#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW(id+1, p, n) - 1 )
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW(id+1, p, n) - BLOCK_LOW(id, p, n))
#define BLOCK_OWNER(index, p, n) (((p) * ((index)+1)-1)/(n))

using namespace std;
int pNum=0, pRank=0;

void processTerminate(double* matrix, double* vector, double* result, double* pResult, double* pMatrixCols, double* pVector, double* pPartialRes, double* pFinal){
	if (pRank == 0){delete [] matrix;delete [] vector;delete [] result;}
	delete [] pMatrixCols;
	delete [] pResult;
	delete [] pVector;
	delete [] pPartialRes;
	delete [] pFinal;
}

void inputInit(int &rows, int &cols, int &size){
	cout << "Enter the size of the matrix and vector in the form of rows, cols\n";
	cout << "Number of rows in Matrix=";
	cin >> rows;
	cout << "Number of columns in the Matrix=";
	cin >> cols;
	cout << "Number of elements in Vector=";
	cin >> size;
	if (size != cols) {
		cout << "The input vector is incompatible with the matrix. Please correct input\n";
		exit(1);
	}
	if (size < pNum) {
		cout << "The size of vector is less than number of processes\n";
	}
}

void randomInit(double* &matrix, double* &vector, int rowCount, int colCount, int size){
	unsigned int seed;
	for (int i = 0; i < rowCount; i++){
		for (int j=0; j < colCount; j++){
			matrix[i * colCount + j] = double(rand_r(&seed) % 100);
		}
	}
	for (int i = 0; i < size; i++){
		vector[i] = double(rand_r(&seed) % 100);
	}
}

void printMatrix(double* matrix, int rowCount, int colCount) {
	for (int i = 0; i < rowCount; i++){
		for (int j = 0; j < colCount; j++){
			cout << matrix[i*colCount + j] << "\t";
		}
		cout << "\n";
	}
	cout << endl;
}

void printVector(double* vector, int size) {
	for (int i = 0; i < size; i++){
		cout << vector[i] << "\t";
	}
	cout << endl;
}

void processInit(double* &matrix, double* &vector, double* &pVector, double* &result, double* &pResult, double* &pPartialRes, double* &pFinal, double* &pMatrixCols, int &rows, int &cols, int &size, int &pNumCols, int &pNumRows) {
	MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	pNumCols = BLOCK_SIZE(pRank, pNum, cols);
	pNumRows = BLOCK_SIZE(pRank, pNum, rows);
	//Memory Allocation in each process
	pVector = new double[pNumCols];
	pMatrixCols = new double [rows * pNumCols];
	pResult = new double[rows];
	pPartialRes = new double [pNumRows * pNum];
	pFinal = new double [pNumRows];

	//Setting the values in the matrix and the vector
	// Input only taken by Proc 0 and then distributed
	if (pRank == 0) {
		vector = new double [size];
		matrix = new double [rows*cols];
		result = new double [rows];
		randomInit(matrix, vector, rows, cols, size);
	}
}

void dataDistribution(double* matrix, double* vector, double* pVector,double* pMatrixCols, int rows, int cols, int size, int pNumCols) {
	int* pSendNum = new int [pNum];
	int* pSendInd = new int [pNum];
	double* row = new double [cols];

	pSendNum[0] = BLOCK_SIZE(0, pNum, cols);
	pSendInd[0] = 0;
	for (int i = 1; i < pNum; i++){
		pSendInd[i] = pSendInd[i-1] + pSendNum[i-1];
		pSendNum[i] = BLOCK_SIZE(i, pNum, cols);
	}
	double* newRow = new double [pSendNum[pRank]];
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			if (pRank == 0) row[j] = matrix[i*cols + j];
		}
		MPI_Scatterv(row, pSendNum, pSendInd, MPI_DOUBLE, newRow, pSendNum[pRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

		for (int j = 0; j < pSendNum[pRank]; j++) {
			pMatrixCols[i*pSendNum[pRank]+j] = newRow[j];
		}
	}
	MPI_Scatterv(vector, pSendNum, pSendInd, MPI_DOUBLE, pVector, pSendNum[pRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	delete [] pSendNum;
	delete [] pSendInd;
	delete [] row;
	delete [] newRow;
}

void colStripeMatrixVectorMul(double* pMatrixCols,double* pVector,int rows,int pNumCols,double* pResult){
	for (int i = 0; i < rows; i++) {
		pResult[i] = 0;
		for (int j = 0; j < pNumCols; j++) {
			pResult[i] += pMatrixCols[i*pNumCols+j] * pVector[j];
		}
	}
}

void partialResultExchange(double* &pFinal, double* pResult, double* pPartialRes, int rows, int pNumRows) {
	int* pSendNum = new int [pNum];
	int* pSendInd = new int [pNum];
	int* pRecvNum = new int [pNum];
	int* pRecvInd = new int [pNum];
	pSendNum[0] = BLOCK_SIZE(0, pNum, rows);
	pSendInd[0] = 0;
	pRecvNum[0] = BLOCK_SIZE(pRank, pNum, rows);
	pRecvInd[0] = 0;
	for (int i = 1; i < pNum; i++){
		pSendInd[i] = pSendInd[i-1] + pSendNum[i-1];
		pSendNum[i] = BLOCK_SIZE(i, pNum, rows);
		pRecvInd[i] = pRecvInd[i-1] + pRecvNum[i-1];
		pRecvNum[i] = BLOCK_SIZE(pRank, pNum, rows);
	}
	MPI_Alltoallv(pResult, pSendNum, pSendInd, MPI_DOUBLE, pPartialRes, pRecvNum, pRecvInd, MPI_DOUBLE, MPI_COMM_WORLD);
	for (int j = 0; j < pNumRows; j++) {
		pFinal[j] = 0.0;
		for (int i = 0; i < pNumRows*pNum; i=i+pNumRows) {
			pFinal[j] += pPartialRes[i+j];
		}
	}
	delete [] pSendNum;
	delete [] pSendInd;
	delete [] pRecvNum;
	delete [] pRecvInd;
}
void resultGather(double* result, double* pFinal, int pNumRows, int rows) {
	int* pReceiveNum = new int [pNum];
	int* pReceiveInd = new int [pNum];
	pReceiveNum[0] = BLOCK_SIZE(0, pNum, rows);
	pReceiveInd[0] = 0;
	for (int i = 1; i < pNum; i++) {
		pReceiveNum[i] = BLOCK_SIZE(i, pNum, rows);
		pReceiveInd[i] = pReceiveInd[i-1] + pReceiveNum[i-1];
	}
	MPI_Gatherv(pFinal, pReceiveNum[pRank], MPI_DOUBLE, result, pReceiveNum, pReceiveInd, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	//MPI_Allgatherv(pFinal, pReceiveNum[pRank], MPI_DOUBLE, result, pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
	delete [] pReceiveNum;
	delete [] pReceiveInd;
}
void matrixVectorMulSerial(double* result, double* matrix, double* vector, int rowCount, int colCount) {
	for (int i = 0; i < rowCount; i++) {
		result[i] = 0;
		for (int j = 0; j < colCount; j++) {
			result[i] += matrix[i*colCount + j]*vector[j];
		}
	}
}

void testResult(double* matrix, double* vector, double* result, int rows, int cols, double &diff_serial) {
	int equal = 0;
	double* serialResult = new double [rows];
	struct timespec start_serial, end_serial;
	clock_gettime(CLOCK_MONOTONIC_RAW, &start_serial);
	matrixVectorMulSerial(serialResult, matrix, vector, rows, cols);
	clock_gettime(CLOCK_MONOTONIC_RAW, &end_serial);
	diff_serial = double((end_serial.tv_sec - start_serial.tv_sec) + (end_serial.tv_nsec - start_serial.tv_nsec)/1000000000.0);
	for (int i = 0; i < rows; i++) {
		if (result[i] != serialResult[i]) equal = 1;
	}
	if (equal == 1) printf("Results of parallel and Serial algorithm are not identical. Please check code\n");
	else printf("The result of Serial and Parallel Algorithm are the same. Good Job! \n ");
	delete [] serialResult;
}
int main(int argc, char* argv[]){
	int rows, cols, size, pNumCols, pNumRows;
	double *matrix,*vector,*result, *pResult, *pMatrixCols, *pVector, *pPartialRes, *pFinal;
 	double start, end;
	double diff_parallel=0.0, diff_serial=0.0, speedup=0.0, average = 0.0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &pNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &pRank);

	rows = stoi(argv[1]); cols = stoi(argv[2]); size = stoi(argv[3]);
	if (pRank == 0) {
		printf("Column-wise matrix vector multiplication\n");
		if ((size != cols) || (size < pNum)) {
			cout << "The input vector is incompatible with the matrix. Or The size of vector is less than number of processes. Please correct input\n";
			exit(1);
		}
	}
	processInit(matrix, vector, pVector, result, pResult,pPartialRes, pFinal, pMatrixCols, rows, cols, size, pNumCols, pNumRows);
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i = 0; i < 100; i++){
		diff_parallel =0.0;
		start = MPI_Wtime();
		dataDistribution(matrix, vector, pVector, pMatrixCols, rows, cols, size, pNumCols);
		colStripeMatrixVectorMul(pMatrixCols, pVector, rows, pNumCols, pResult);
		partialResultExchange(pFinal, pResult, pPartialRes, rows, pNumRows);
		resultGather(result, pFinal, pNumCols, rows);
		end = MPI_Wtime();
		diff_parallel = end - start;
		average += diff_parallel;
	}
	average /= 100;
	if (pRank == 0) testResult(matrix, vector, result, rows, cols, diff_serial);
	if (pRank == 0) {
		speedup = diff_serial/average;
		printf("rows = %d\ncols = %d\nSerial Time = %fs\nMPI Time = %fs\nSpeedup = %f\n", rows, cols, diff_serial, average, speedup);
	}
	processTerminate(matrix, vector, result, pResult, pMatrixCols, pVector, pPartialRes, pFinal);
	MPI_Finalize();
  return 0;
}
