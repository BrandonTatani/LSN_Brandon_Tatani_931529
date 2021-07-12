#include <iostream>
#include "random.h"
#include "Genetic_Algorithm.h"
using namespace std;

int main(){

//  Random Generator
	Random* RNG = new Random(1);
	int n_city = 32;
	double raggio = 10.;
	int popolazione = 500;
	int generazioni = 200;
	int iprint = generazioni/10, nprint=0;
	City* circo = new City(n_city, raggio, RNG);
	
	Genetic root(popolazione, circo, RNG);
	
	circo -> Print();
	
	for(int i=0; i<generazioni; i++){

		root.Update_Param();
		root.Ordine();

		root.Print_Best();
		root.half_ave();
		
		root.Crossover();
		root.Mutazione();

		if(i%iprint == 0)	cout << nprint++ << "0%" << endl;
		
	}
	
	root.Update_Param();
	root.Ordine();
	root.Print_Cromosoma();
	
	cout << "100%" << endl;
	
	RNG -> SaveSeed();
	
return 0;
}
