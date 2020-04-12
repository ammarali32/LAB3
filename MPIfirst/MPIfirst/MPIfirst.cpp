#include<mpi.h>
#include <iostream>
#include <math.h>
#include <chrono>
#include <fstream>
#include <string>

const int MAX_SIZE = 500000;
int size, myRank;

int PARTITION(int* A, int l, int r) {
	int pi = A[l];
	int i = l - 1, j = r + 1;
	while (1) {
		do {
			i++;
		}while (A[i] < pi);
		do {
			j--;
		}while (A[j] > pi);
		if (i >= j)return j;
		int temp = A[i];
		A[i] = A[j];
		A[j] = temp;
	}
}

int Average(int* A, int l, int r) {
	long long sum = 0;
	for (int i = l; i <= r; i++) {
		sum += A[i];
	}
	return (int)(sum / (r - l + 1));
}
int LOM_PARTITION(int* A, int l, int r, int pi) {
	int i = l - 1;
	for (int j = l; j < r; j++) {
		if (A[j] <= pi) {
			i++;
			int temp = A[i];
			A[i] = A[j];
			A[j] = temp;
		}
	}
	int temp = A[i + 1];
	A[i + 1] = A[r];
	A[r] = temp;
	return (A[i + 1] > pi) ? i : (i + 1);
}

void quicksort(int* A, int l, int r) {
	if (l < r) {
		int pi = PARTITION(A, l, r);
		quicksort(A, l, pi);
		quicksort(A, pi + 1, r);
	}
}

int* Init(int& nPar) {
	if (myRank == 0) {
		int* Input = (int*)malloc(MAX_SIZE * sizeof(int));
		freopen("input.txt", "r", stdin);
		int n;
		std::cin >> n;
		for (int i = 0; i < n; i++) {
			std::cin >> Input[i];
		}
		nPar = n / size;
		int* partData = (int*)malloc(nPar * sizeof(int));
		for (int i = 0; i < nPar; i++) {
			partData[i] = Input[i];
		}
		for (int rank = 1; rank < size; rank++) {
			MPI_Send(&nPar, 1, MPI_INT, rank, 0, MPI_COMM_WORLD);
			MPI_Send((Input + rank * nPar), nPar, MPI_INT, rank, 1, MPI_COMM_WORLD);
		}
		free(Input);
		return partData;
	}
	else {
		MPI_Recv(&nPar, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		int* partData = (int*)malloc(nPar * sizeof(int));
		MPI_Recv(partData, nPar, MPI_INT, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		return partData;
	}
}

void FindPivot(int* piPartData, int nPar, int dimension, int mask, int it, int* piarray) {
	MPI_Group WorldGroup;
	MPI_Group cubeGroup;
	MPI_Comm cubeComm;
	int groupNum = size / (int)pow(2, dimension - it);
	int* pRanks = (int*)malloc(groupNum * sizeof(int));
	int StartProc = myRank - groupNum;
	StartProc = ((StartProc < 0) ? 0 : StartProc);
	int EndProc = myRank + groupNum;
	EndProc = ((EndProc > size) ? size : EndProc);
	int j = 0;
	
	for (int p = StartProc; p < EndProc; p++) {
		//std::cout << myRank << " " << mask << " " << p << " " << it << std::endl;
		if ((myRank & mask) >> (it) == (p & mask) >> (it)) {
			pRanks[j++] = p;
			//std::cout << pRanks[j - 1] << std::endl;
		}
		//pRanks[j++] = p;
	}
	
	MPI_Comm_group(MPI_COMM_WORLD, &WorldGroup);
	MPI_Group_incl(WorldGroup, groupNum, pRanks, &cubeGroup);
	MPI_Comm_create(MPI_COMM_WORLD, cubeGroup, &cubeComm);
	if (myRank == pRanks[0]) {
		*piarray = Average(piPartData, 0, nPar - 1);
	}
	
	MPI_Bcast(piarray, 1, MPI_INT, 0, cubeComm);
	MPI_Group_free(&cubeGroup);
	MPI_Comm_free(&cubeComm);
	free(pRanks);
}
void MERGE(int* mergedData, int mergDataSize, int* partData, int DataSize, int * partRecData, int RecDataSize) {
	for (int i = 0; i < mergDataSize; i++) {
		mergedData[i] =( (i < DataSize) ? partData[i] : partRecData[i - DataSize]);
	}
}

void PARALLELQUICKSORT(int *& partData,int * nPar){
	MPI_Status status;
	int CommRank;
	int* PtData, * partSendData;
	int DataSize, SendingDataSize, RecDataSize, MergDataSize;
	int CubeDim = (int)(log(size) / log(2));
	int mask = size;
	int pi;
	for (int i = CubeDim; i > 0; i--) {
		FindPivot(partData, *nPar, CubeDim, mask, i, &pi);
		mask = mask >> 1;
		int id = LOM_PARTITION(partData, 0, *nPar - 1, pi);
		if ((myRank & mask) == 0) {
			CommRank = myRank + mask;
			partSendData = &partData[id + 1];
			SendingDataSize = *nPar - id - 1;
			if (SendingDataSize < 0)SendingDataSize = 0;
			PtData = &partData[0];
			DataSize = id + 1;
		}
		else {
			CommRank = myRank - mask;
			partSendData = &partData[0];
			SendingDataSize = id + 1;
			if (SendingDataSize > *nPar)SendingDataSize = id;
			PtData = &partData[id + 1];
			DataSize = *nPar - id - 1;
			if (DataSize < 0) DataSize = 0;
		}
		MPI_Sendrecv(&SendingDataSize, 1, MPI_INT, CommRank, 0, &RecDataSize, 1, MPI_INT, CommRank, 0, MPI_COMM_WORLD, &status);
		int* pRecvData = (int*)malloc(RecDataSize * sizeof(int));
		MPI_Sendrecv(partSendData, SendingDataSize, MPI_INT, CommRank, 0, pRecvData, RecDataSize, MPI_INT, CommRank, 0, MPI_COMM_WORLD, &status);
		MergDataSize = DataSize + RecDataSize;
		int* pMergeData = (int*)malloc(MergDataSize * sizeof(int));
		MERGE(pMergeData, MergDataSize, PtData, DataSize, pRecvData, RecDataSize);
		partData = pMergeData;
		*nPar = MergDataSize;
		free(pRecvData);
	}
}
void ProcessTermination(int* ProcData, int nPar) {

	MPI_Status status;
	long long DataSize = 0;
	MPI_Reduce(&nPar, &DataSize, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	int* Data, * rcounts, * displs;
	rcounts = (int*)malloc(size * sizeof(int));
	Data = (int*)malloc(DataSize * sizeof(int));
	displs = (int*)malloc(size * sizeof(int));
	MPI_Gather(&nPar, 1, MPI_INT, rcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (myRank == 0) {
		int acc = 0;
		for (int r = 0; r < size; r++) {
			displs[r] = acc;
			acc += rcounts[r];
		}
	}
	MPI_Gatherv(ProcData, nPar, MPI_INT, Data, rcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
	if (myRank == 0)
	{
		freopen("output.txt", "w", stdout);
		std::cout << DataSize << std::endl;
		for (int i = 0; i < DataSize; i++) {
			std::cout << Data[i] << " ";
		}
		
	}
	free(ProcData);
	
}


int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int ProcDataSize;
	int* ProcData = Init( ProcDataSize);
	auto start = std::chrono::high_resolution_clock::now();
	PARALLELQUICKSORT(ProcData, &ProcDataSize);
	quicksort(ProcData, 0, ProcDataSize - 1);
	auto finish = std::chrono::high_resolution_clock::now();
	
	ProcessTermination( ProcData, ProcDataSize);
	std::cout << (std::chrono::duration<double, std::milli>(finish - start).count()) << std::endl;
	MPI_Finalize();
	return 0;
}

