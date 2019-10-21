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
int pNum=0, pRank=0, gridId;

void processTerminate(double* matrix, double* vector, double* result, double* pResult, double* pMatrix, double* pVector) {
	if (pRank == 0){delete [] matrix; delete [] vector; delete [] result;}
	delete [] pMatrix;
	delete [] pResult;
	delete [] pVector;
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
	cout.flush();
}
void printVector(double* vector, int size) {
	for (int i = 0; i < size; i++){
		cout << vector[i] << "\t";
	}
	cout << endl;
	cout.flush();
}
void processInit(double* &matrix, double* &vector, double* &pMatrix, double* &pVector, double* &result, double* &pResult, int rows, int cols, int size, int &pNumCols, int &pNumRows, int gridDims, int* gridSize, int* gridCoord, MPI_Comm checkerboard_comm) {
	MPI_Bcast(&rows, 1, MPI_INT, 0, checkerboard_comm);
	MPI_Bcast(&cols, 1, MPI_INT, 0, checkerboard_comm);
	MPI_Bcast(&size, 1, MPI_INT, 0, checkerboard_comm);
	int gridPeriod[gridDims];
	MPI_Cart_get(checkerboard_comm, gridDims, gridSize, gridPeriod, gridCoord);
	pNumRows = BLOCK_SIZE(gridCoord[0], gridSize[0], rows);
	pNumCols = BLOCK_SIZE(gridCoord[1], gridSize[1], cols);
	pMatrix = new double [pNumRows*pNumCols];
	pVector = new double [pNumCols];
	pResult = new double[pNumRows];

	//Setting the values in the matrix and the vector
	// Input only taken by Proc 0 and then distributed
	if (gridId == 0) {
		vector = new double [size];
		matrix = new double [rows*cols];
		result = new double [rows];
		randomInit(matrix, vector, rows, cols, size);
	}
}
void matrixDistribution(double* matrix, double* pMatrix, int rows, int cols , int pNumRows, int pNumCols, int gridDims, int* gridSize, int* gridCoord, MPI_Comm checkerboard_comm) {
	int coords[gridDims], localRows, localCols, dest_id;
	MPI_Status status;
	double* oneRow = new double [cols];
	for (int i = 0; i < gridSize[0]; i++) {
		coords[0] = i;
		localRows = BLOCK_SIZE(i, gridSize[0], rows);
		for (int j = 0; j < localRows; j++) {
			if (gridId == 0){
				for (int r = 0; r < cols; r++) {
					oneRow[r] = matrix[(i*localRows + j)*cols + r];
				}
			}
			for (int k = 0; k < gridSize[1]; k++){
				coords[1] = k;
				localCols = BLOCK_SIZE(k, gridSize[1], cols);
				MPI_Cart_rank(checkerboard_comm, coords, &dest_id);
				if (gridId == 0){
					if (dest_id == 0){
						for (int r = 0; r < localCols; r++){
							pMatrix[j*localCols+r] = oneRow[r];
						}
					}else {
						void* sendRow;
						sendRow = oneRow + (k*localCols);
						MPI_Send(sendRow, localCols, MPI_DOUBLE, dest_id, 0, checkerboard_comm);
					}
				}else if (gridId == dest_id){
					double* newRow = new double [localCols];
					MPI_Recv(newRow, localCols, MPI_DOUBLE, 0, 0, checkerboard_comm, &status);
					for (int r = 0; r < localCols; r++){
						pMatrix[j*localCols+r] = newRow[r];
					}
					delete [] newRow;
				}
			}
		}
	}
	delete [] oneRow;
}
void vectorDistribution(double* vector, double* pVector, int rows, int cols, int pNumCols, int gridDims, int* gridSize, int* gridCoord, MPI_Comm checkerboard_comm, MPI_Comm row_comm, MPI_Comm col_comm){
	int coords_src[gridDims], coords_dest[gridDims], src_id, dest_id;
	MPI_Status status;

	if (gridSize[0] == gridSize[1]){
		if (gridCoord[1] == 0){
			int* pSendNum = new int [gridSize[0]];
			int* pSendInd = new int [gridSize[0]];
			pSendNum[0] = BLOCK_SIZE(0, gridSize[0], cols);
			pSendInd[0] = 0;
			for (int i = 0; i < gridSize[0]; i++) {
				pSendInd[i] = pSendInd[i-1] + pSendNum[i-1];
				pSendNum[i] = BLOCK_SIZE(i, gridSize[0], cols);
			}
			int cRank;
			MPI_Comm_rank(col_comm, &cRank);
			double* pTempVector = new double [pSendNum[cRank]];
			MPI_Scatterv(vector, pSendNum, pSendInd, MPI_DOUBLE, pTempVector, pSendNum[cRank], MPI_DOUBLE, 0, col_comm);
			coords_dest[0] = gridCoord[1];
			coords_dest[1] = gridCoord[0];
			MPI_Cart_rank(checkerboard_comm, coords_dest, &dest_id);
			if (gridCoord[0] == 0) {
				for (int i =0; i < pNumCols; i++) {
					pVector[i] = pTempVector[i];
				}
			} else MPI_Send(pTempVector, pSendNum[cRank], MPI_DOUBLE, dest_id, 0, checkerboard_comm);
			delete [] pSendNum;
			delete [] pSendInd;
			delete [] pTempVector;
		} else if (gridCoord[0] == 0){
			coords_src[0] = gridCoord[1]; coords_src[1] = gridCoord[0];
			MPI_Cart_rank(checkerboard_comm, coords_src, &src_id);
			MPI_Recv(pVector, pNumCols, MPI_DOUBLE, src_id, 0, checkerboard_comm, &status);
		}
	} else {
		if (gridCoord[0] == 0) {
			int* pSendNum = new int [gridSize[1]];
			int* pSendInd = new int [gridSize[1]];
			pSendNum[0] = BLOCK_SIZE(0, gridSize[1], cols);
			pSendInd[0] = 0;
			for (int i = 0; i < gridSize[1]; i++) {
				pSendInd[i] = pSendInd[i-1] + pSendNum[i-1];
				pSendNum[i] = BLOCK_SIZE(i, gridSize[1], cols);;
			}
			MPI_Scatterv(vector,  pSendNum, pSendInd, MPI_DOUBLE, pVector, pSendNum[gridId], MPI_DOUBLE, 0, row_comm);
			delete [] pSendNum;
			delete [] pSendInd;
		}
	}
	MPI_Bcast(pVector, pNumCols, MPI_DOUBLE, 0, col_comm);
}
void checkerboardMatrixVectorMul(double* pMatrix,double* pVector,int pNumRows,int pNumCols,double* pResult){
	for (int i = 0; i < pNumRows; i++) {
		pResult[i] = 0;
		for (int j = 0; j < pNumCols; j++) {
			pResult[i] += pMatrix[i*pNumCols+j] * pVector[j];
		}
	}
}
void resultGather(double* result, double* pResult, int rows, int pNumRows, int pNumCols, int* gridSize, int* gridCoord, MPI_Comm row_comm, MPI_Comm col_comm) {
	if (gridCoord[1] == 0) {
		MPI_Reduce(MPI_IN_PLACE, pResult, pNumRows, MPI_DOUBLE, MPI_SUM, 0, row_comm);
	} else {
		MPI_Reduce(pResult, pResult, pNumRows, MPI_DOUBLE, MPI_SUM, 0, row_comm);
	}
	if (gridCoord[1] == 0){
		int* pReceiveNum = new int [gridSize[0]];
		int* pReceiveInd = new int [gridSize[0]];
		pReceiveNum[0] = BLOCK_SIZE(0, gridSize[0], rows);
		pReceiveInd[0] = 0;
		for (int i = 1; i < gridSize[0]; i++) {
			pReceiveNum[i] = BLOCK_SIZE(i, gridSize[0], rows);
			pReceiveInd[i] = pReceiveInd[i-1] + pReceiveNum[i-1];
		}
		MPI_Gatherv(pResult, pNumRows, MPI_DOUBLE, result, pReceiveNum, pReceiveInd, MPI_DOUBLE, 0, col_comm);
		delete [] pReceiveInd;
		delete [] pReceiveNum;
	}
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
	double *matrix,*vector,*result, *pResult, *pMatrix, *pVector;
 	double start, end;
	double diff_parallel=0.0, diff_serial=0.0, speedup=0.0, average = 0.0;

	MPI_Init(&argc, &argv);
	rows = stoi(argv[1]); cols = stoi(argv[2]); size = stoi(argv[3]);

	MPI_Comm_size(MPI_COMM_WORLD, &pNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &pRank);

	if (pRank == 0) {
		printf("Checkerboard matrix vector multiplication\n");
		if ((size != cols) || (size < pNum)) {
			cout << "The input vector is incompatible with the matrix. Or The size of vector is less than number of processes. Please correct input\n";
			exit(1);
		}
	}
	//Creating the new communicator for reading checkerboard matrix
	int gridDims=2; int gridSize[gridDims], periodic[gridDims], gridCoord[gridDims];
	MPI_Comm checkerboard_comm;
	MPI_Comm row_comm, col_comm; // Row and Column Commnicators for communicating in one row or column of the grid

	gridSize[0] = 0; gridSize[1] = 0; periodic[0] = periodic[1] =0; // No wrap-around in the processor grid
	MPI_Dims_create(pNum, gridDims, gridSize); //returns grid size like 4x4 using dimensions.
  MPI_Cart_create(MPI_COMM_WORLD, gridDims, gridSize, periodic, 1, &checkerboard_comm);
	MPI_Comm_rank(checkerboard_comm, &gridId);
	MPI_Comm_size(checkerboard_comm, &pNum);
	processInit(matrix, vector, pMatrix, pVector, result, pResult, rows, cols, size, pNumCols, pNumRows, gridDims, gridSize, gridCoord, checkerboard_comm);
	MPI_Comm_split(checkerboard_comm, gridCoord[0], gridCoord[1], &row_comm);
	MPI_Comm_split(checkerboard_comm, gridCoord[1], gridCoord[0], &col_comm);
	for (int i = 0; i < 100; i++) {
		diff_parallel = 0.0;
		start = MPI_Wtime();
		matrixDistribution(matrix,pMatrix, rows, cols, pNumRows, pNumCols, gridDims, gridSize, gridCoord, checkerboard_comm);
		vectorDistribution(vector, pVector, rows, cols, pNumCols, gridDims, gridSize, gridCoord, checkerboard_comm, row_comm, col_comm);
		checkerboardMatrixVectorMul(pMatrix, pVector, pNumRows, pNumCols, pResult);
		resultGather(result, pResult, rows, pNumRows, pNumCols, gridSize, gridCoord, row_comm, col_comm);
		end = MPI_Wtime();
		diff_parallel = end - start;
		average += diff_parallel;
	}
	average /= 100;
	if (gridId == 0) testResult(matrix, vector, result, rows, cols, diff_serial);
	if (gridId == 0) {
		speedup = diff_serial/average;
		printf("rows = %d\ncols = %d\nSerial Time = %fs\nMPI Time = %fs\nSpeedup = %f\n", rows, cols, diff_serial, average, speedup);
	}

	processTerminate(matrix, vector, result, pResult, pMatrix, pVector);
	MPI_Finalize();
  return 0;
}
