#include<iostream>
#include<math.h>
#include<chrono>
#include<fstream>
#include"mpi.h"
#include<stdio.h>
#include<stdlib.h>
#include<algorithm>
#define MAX 1000005
using namespace std;
int LOM_PARTITION(int* A, int l, int r) {
	int pi = A[r];
	int j = l - 1;
	for (int i = l; i < r; i++) {
		if (A[i] <= pi) {
			j++;
			std::swap(A[i], A[j]);
		}
	}
	std::swap(A[j + 1], A[r]);
	return j + 1;
}

void NORM_QUICK_SORT(int* A, int l, int r) {
	
	if (l < r) {
		int par = LOM_PARTITION(A, l, r);
		NORM_QUICK_SORT(A, l, par-1);
		NORM_QUICK_SORT(A, par + 1, r);
	}
}
int SORT_REC(int* A, int Asize, int process, int maxprocesses, int depth) {
	MPI_Status status;
	int sharing = process + (1 << depth);
	depth++;
	//std::cout << sharing <<" "<<process<<" "<<depth<<" "<<maxprocesses<< std::endl;
	if (sharing > maxprocesses) {
		NORM_QUICK_SORT(A, 0, Asize-1 );
		return 0;
	}
	int id = 0, piv_id;
	do {
		piv_id = LOM_PARTITION(A, id, Asize - 1);
		id++;
	} while (piv_id == id - 1 && id <= Asize -1);
	if (id > Asize - 1) {
		return 0;
	}
	if (piv_id <= Asize - piv_id) {
		MPI_Send(A, piv_id - 1, MPI_INT, sharing, piv_id, MPI_COMM_WORLD);
		SORT_REC((A + piv_id + 1), (Asize - piv_id - 1), process, maxprocesses, depth);
		MPI_Recv(A, piv_id - 1, MPI_INT, sharing, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	}
	else {
		MPI_Send((A + piv_id + 1), Asize - piv_id - 1, MPI_INT, sharing, piv_id + 1, MPI_COMM_WORLD);
		SORT_REC(A, piv_id + 1, process, maxprocesses, depth);
		MPI_Recv((A + piv_id + 1), Asize - piv_id - 1, MPI_INT, sharing, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
	}
}
int A[MAX],B[MAX];

int main(int argc,char *argv[]){
	ifstream inf("test.csv");

	char c;
	string ss; ss.clear();
	while (inf.get(c))
	{
		ss.push_back(c);
	}
	inf.close();
	ofstream of("test.csv");
	of << ss << endl; 
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int depth = 0;
	while ((1 << depth) <= rank)depth++;
	int n = 10;
	of << " using " << size << " of threads\n";
	while (n <= 1000000) {
		if (rank == 0) {
			
			/*freopen("input.txt", "r", stdin);

			int n; std::cin >> n;
			for (int i = 0; i < n; i++)std::cin >> A[i];*/
			//int n = 100000;
			//for (int i = 0; i < n; i++)A[i] = rand() % 100 +1;

			for (int i = 0; i < n; i++) {
				A[i] = rand() % 10000 + 1;
			}
			cout << n << " stopped here " << rank << " " << size << endl;
			auto start = std::chrono::steady_clock::now();
			SORT_REC(A, n, rank, size - 1, depth);
			auto finish = std::chrono::steady_clock::now();
			of << (std::chrono::duration<double, std::milli>(finish - start).count()) << " ms " << std::endl;
			n *= 10;

			//freopen("output.txt", "w", stdout);
			//std::cout << n << std::endl;
			//for (int i = 0; i < n; i++)std::cout << A[i]<<" ";
			//fclose(stdout);			
		}
		else {
			MPI_Status status;
			int arsize;
			MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			MPI_Get_count(&status, MPI_INT, &arsize);
			int source = status.MPI_SOURCE;
			int* ar = new int[arsize];
			MPI_Recv(ar, arsize, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			SORT_REC(ar, arsize, rank, size - 1, depth);
			MPI_Send(ar, arsize, MPI_INT, source, 0, MPI_COMM_WORLD);
		}
		
		
	}
	of.close();
	MPI_Finalize();
	return 0;

	
	
}


