#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

// Creazione generatore di Random	
  
	Random rnd(0);

// Variabili dell'applicazione

	int M = 100000, N=100, L = M/N;
	vector <double> f(N), sigma2(N), n(N);
	double sum = 0;
	double sum_prog, sum2_prog;
	

//	Calcolo di f[i] per ogni blocco
				
	for(int i=0; i<N; i++){
		sum = 0;
		
		for(int j=0; j<L; j++)
			sum += Coseno(rnd.Rannyu());
			
		f[i] = sum/L;
	}
	
//	Calcolo medie e incertezze progressive di I

	ofstream file("../output2_11.csv");

	file << "I" << "," << "errore" << endl;

	for(int i=0; i<N; i++){
		sum_prog = 0;
		sum2_prog = 0;

		for(int j=0; j<=i; j++){
			sum_prog += f[j];
			sum2_prog += pow(f[j], 2);
		}
		sum_prog /= (i+1);
		sum2_prog /= (i+1);
		
		file << sum_prog-1. << "," << error(sum_prog, sum2_prog, i) << endl;
		
	} 
	
	file.close();	
	

//	Calcolo di f[i] per ogni blocco
				
	for(int i=0; i<N; i++){
		sum = 0;
		
		for(int j=0; j<L; j++)
			sum += Coseno(rnd.Rannyu());
			
		f[i] = sum/L;
	}
	
/////////////////////// -> Con distribuzione affine <- /////////////////////////////////
	
//	Calcolo di f[i] per ogni blocco
	double tmp;
				
	for(int i=0; i<N; i++){
		sum = 0;
		
		for(int j=0; j<L; j++){
			tmp = rnd.Reject();
			sum += 5*M_PI/12. * cos(M_PI*0.5*tmp) / (1-tmp*tmp*0.5);
		}
		f[i] = sum/L;
	}



//	Calcolo medie e incertezze progressive di I

	file.open("../output2_12.csv");

	file << "I" << "," << "errore" << endl;

	for(int i=0; i<N; i++){
		sum_prog = 0;
		sum2_prog = 0;

		for(int j=0; j<=i; j++){
			sum_prog += f[j];
			sum2_prog += pow(f[j], 2);
		}
		sum_prog /= (i+1);
		sum2_prog /= (i+1);
		
		file << sum_prog-1. << "," << error(sum_prog, sum2_prog, i) << endl;
		
	} 
	
	file.close();	



// Salvataggio del seme di uscita
   rnd.SaveSeed();
   return 0;
}

