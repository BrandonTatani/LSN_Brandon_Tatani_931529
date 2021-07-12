#include <iostream>
#include "random.h"
#include "Simulated_Annealing.h"
using namespace std;

int main(){

//  Random Generator
	Random* RNG = new Random(1);
	int n_city = 32;
	double raggio = 10.;
	double T_in = 2.5, T_fin = 0.1;
	int step = 2000;
	double h = 0.1;
	City* circo = new City(n_city, raggio, RNG);
	
	Annealing root(circo, RNG, T_in);
	
	circo -> Print();
	
	for(double i=T_in; i>T_fin; i=i-h){
		root.SetT(i);
		
		for(int k=0; k<step; k++){
			root.Step();
			root.Save_Param();

		}
		
		cout << i << endl;
	}
	
	root.Update_Param();
	root.Save_Percorso();
	RNG -> SaveSeed();
	
return 0;
}
