 #include <iostream>
#include "random.h"
#include "mpi.h"
#include "Genetic_Algorithm.h"
#include <vector>
#include <cmath>
using namespace std;

int main(int argc, char* argv[]){

//  Random Generator
	
	int n_city = 32;
	double raggio = 10.;
	int popolazione = 1000;
	int generazioni = 1200;
	
	int size, rank, N_migr = 100;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status stat;
	

	
	int iprint = generazioni/10, nprint=0;
	City* circo = new City(n_city);
		
	Random* RNG = new Random(1, rank);
	Genetic root(popolazione, circo, RNG);

	
	/////////////////////////////////////
	
	int* tmp = new int(n_city);
	tmp = root.Get_Address(0,0);
	int itag = 1;
	int imsg = rank;
/*	cout << &tmp[31] << " " << rank << endl;
	
	MPI_Bcast(&tmp[0],n_city, MPI_INTEGER,1,MPI_COMM_WORLD);
	
	if(rank == 1){
		MPI_Send(&tmp[0],n_city, MPI_INTEGER,0,itag, MPI_COMM_WORLD);
	}else if(rank ==  0){
		MPI_Recv(&tmp[0],n_city, MPI_INTEGER,1,itag,MPI_COMM_WORLD, &stat);
	}
//	cout << imsg << endl;
	cout << &tmp[31] << " " << rank << endl;
	
*/	
	/////////////////////////////////////
	
	/*cout << root.Get_Address(0, 0) << endl;
	
	int* test = new int(32);
	test = root.Get_Address(0,0);
	cout << test[0] << " " << root.Get_Gene(0,0) << endl;
	*/
		
//	circo -> Print();
	
	for(int i=0; i<generazioni; i++){
		root.Update_Param();
		root.Ordine();
		/////////////////////
		if(i%N_migr==0){
			for(int k=0; k<size; k++){
				itag = int(popolazione * pow(RNG->Rannyu(), 2.5));
				tmp = root.Get_Address(0,0);
				MPI_Bcast(&tmp[0],n_city, MPI_INTEGER,k,MPI_COMM_WORLD);
				cout << rank << " " << root.GetL(itag) << endl;
				if(rank != k)
					root.Copy(tmp, popolazione-1-k);	
			}
//			cout << "Migrazione " << rank << endl;
//			if(rank == 0) root.Print(popolazione-2);
//			if(rank == 1) root.Print(0);
			root.Update_Param();
			root.Ordine();
		}
		/////////////////////
		
		if(rank == 0){
			root.Print_Best();
			root.half_ave();
		}
		
		root.Crossover();
		root.Mutazione();
		
		if(rank == 0 && i%iprint == 0)	cout << nprint++ << "0%" << endl;
	}
	
	root.Update_Param();
	root.Ordine();
	if(rank==0)	root.Print_Cromosoma();
	
	cout << "100%" << endl;
	if(rank == 0)	RNG -> SaveSeed();
	
	MPI_Finalize();
	
return 0;
}
