#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

// Creazione generatore di Random	
  
	int seed[4] = {0, 0, 2, 1} ;
	Random rnd(seed);

// Variabili dell'applicazione

	int Nhit = 0, M = 1000000, N=100, H = M/N;
	double d = 1., L = 0.5;
	vector <double> pi(N);
	double sum_prog, sum2_prog;
	
	double x, Lcos;  /* Poiché il problema è invariante lungo y
				estraggo casualmente x punto medio
				della sbarra e cos che rappresenta
				la proiezione lungo x dell'ago. */

//	Calcolo di pi per ogni blocco
				
	for(int i=0; i<N; i++){
		Nhit=0;
		for(int j=0; j<H; j++){
			x = rnd.Rannyu(0,10);
			
			Lcos = L/2*rnd.Reject();
		//	Lcos = L/2*Coseno(rnd.Rannyu(0, M_PI));
			if(fabs(Lcos) >= fabs(x-round(x)))
				Nhit++;		
		}
		pi[i] = (2*L/d)*(float(H)/Nhit);
		cout << pi[i] << "\t";
	}
	
//	Calcolo medie e incertezze progressive

	ofstream file("output_test.csv");

	file << "x" << "," << "y" << endl;

	for(int i=0; i<N; i++){
		sum_prog = 0;
		sum2_prog = 0;

		for(int j=0; j<=i; j++){
			sum_prog += pi[j];
			sum2_prog += pow(pi[j], 2);
		}
		sum_prog /= (i+1);
		sum2_prog /= (i+1);
		
		file << sum_prog << "," << error(sum_prog, sum2_prog, i) << endl;
		
	} 
	
	file.close();	


// Salvataggio del seme di uscita
   rnd.SaveSeed();
   return 0;
}

