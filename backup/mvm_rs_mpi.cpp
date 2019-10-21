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

void processTerminate(double* matrix, double* vector, double* result, double* pResult, double* pMatrixRows){
	if (pRank == 0){delete [] matrix;}
	delete [] vector;
	delete [] result;
	delete [] pMatrixRows;
	delete [] pResult;
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
void processInit(double* &matrix, double* &vector, double* &result, double* &pResult, double* &pMatrixRows, int &rows, int &cols, int &size, int &pNumRows) {
	//if (pRank == 0) {
	//	do{
	//	inputInit(rows, cols, size);
	//	}while(size < pNum);
	//}
	MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
	pNumRows = BLOCK_SIZE(pRank, pNum, rows);
	//Memory Allocation in each process
	vector = new double[size];
	result = new double[rows];
	pMatrixRows = new double [pNumRows*cols];
	pResult = new double[pNumRows];

	//Setting the values in the matrix and the vector
	// Input only taken by Proc 0 and then distributed
	if (pRank == 0) {
		matrix = new double [rows*cols];
		randomInit(matrix, vector, rows, cols, size);
		//printMatrix(matrix, rows, cols);
		//printVector(vector, size);
	}
}
void dataDistribution(double* matrix, double* vector, double* pMatrixRows, int rows, int cols, int size, int pNumRows) {
	//Broadcast the vector to all processes
	MPI_Bcast(vector, size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int* pSendNum = new int [pNum];
	int* pSendInd = new int [pNum];
	pSendNum[0] = BLOCK_SIZE(0, pNum, rows) * cols;
	pSendInd[0] = 0;
	for (int i = 1; i < pNum; i++){
		pNumRows = BLOCK_SIZE(i, pNum, rows);
		pSendInd[i] = pSendInd[i-1] + pSendNum[i-1];
		pSendNum[i] = pNumRows * cols;
	}
	if (pRank == 0) {
		for (int i = 0; i < pNum; i++){
			//cout << pSendNum[i] << "\t" << pSendInd[i] << endl;
		}
	}
	MPI_Scatterv(matrix, pSendNum, pSendInd, MPI_DOUBLE, pMatrixRows, pSendNum[pRank], MPI_DOUBLE, 0, MPI_COMM_WORLD);
	delete [] pSendNum;
	delete [] pSendInd;
}
void rowStripeMatrixVectorMul(double* pMatrixRows,double* vector,int pNumRows,int cols,double* pResult){
	for (int i = 0; i < pNumRows; i++) {
		pResult[i] = 0;
		for (int j = 0; j < cols; j++) {
			pResult[i] += pMatrixRows[i*cols+j] * vector[j];
		}
	}
}
void resultGather(double* result, double* pResult, int pNumRows, int rows) {
	int* pReceiveNum = new int [pNum];
	int* pReceiveInd = new int [pNum];
	pReceiveNum[0] = BLOCK_SIZE(0, pNum, rows);
	pReceiveInd[0] = 0;
	for (int i = 1; i < pNum; i++) {
		pReceiveNum[i] = BLOCK_SIZE(i, pNum, rows);
		pReceiveInd[i] = pReceiveInd[i-1] + pReceiveNum[i-1];
	}
	if (pRank == 0) {
		for (int i = 0; i < pNum; i++){
			//cout << pReceiveNum[i] << "\t" << pReceiveInd[i] << endl;
		}
	}
	MPI_Allgatherv(pResult, pReceiveNum[pRank], MPI_DOUBLE, result, pReceiveNum, pReceiveInd, MPI_DOUBLE, MPI_COMM_WORLD);
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
void testResult(double* matrix, double* vector, double* result, int rows, int cols) {
	if (pRank == 0) {
		int equal = 0;
		double* serialResult = new double [rows];
		matrixVectorMulSerial(serialResult, matrix, vector, rows, cols);
		for (int i = 0; i < rows; i++) {
			if (result[i] != serialResult[i]) equal = 1;
		}
		if (equal == 1) printf("Results of parallel and Serial algorithm are not identical. Please check code\n");
		else printf("The result of Serial and Parallel Algorithm are the same. Good Job! \n ");
	}

}
int main(int argc, char* argv[]){
	int rows, cols, size, pNumRows;
	double *matrix,*vector,*result, *pResult, *pMatrixRows;
 	double start, end;
	double diff_parallel=0.0, diff_serial=0.0, speedup=0.0;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &pNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &pRank);

	rows = stoi(argv[1]); cols = stoi(argv[2]); size = stoi(argv[3]);
	if ((size != cols) || (size < pNum)) {
		cout << "The input vector is incompatible with the matrix. Or The size of vector is less than number of processes. Please correct input\n";
		exit(1);
	}
	if (pRank == 0) {
		printf("Row-wise matrix vector multiplication\n");
	}
	processInit(matrix, vector, result, pResult, pMatrixRows, rows, cols, size, pNumRows);
	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();
	dataDistribution(matrix, vector, pMatrixRows, rows, cols, size, pNumRows);
	rowStripeMatrixVectorMul(pMatrixRows, vector, pNumRows, cols, pResult);
	resultGather(result, pResult, pNumRows, rows);
	end = MPI_Wtime();
	diff_parallel = end - start;
	testResult(matrix, vector, result, rows, cols);
	//if (pRank == 7) { cout << "Rows = \t" << pNumRows << endl;}
	if (pRank == 0) {
		printf ("Duration of parallel algorithm =%f\n", diff_parallel);
		//printMatrix(matrix, rows, cols);
	}
	if (pRank == 6) {
		//printMatrix(pMatrixRows, pNumRows, cols);
		//printVector(result, rows);
	}
	//cout << pNum << "\t" << pRank << "\n";
	MPI_Finalize();
//	processTerminate(matrix, vector, result, pResult, pMatrixRows);
  return 0;
}
