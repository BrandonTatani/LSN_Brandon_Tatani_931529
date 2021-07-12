#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

// Creazione generatore di Random	
  
	Random rnd(1);

// Variabili dell'applicazione

	int M = 100000, N=100, L = M/N;
	vector <double> r(N), sigma2(N);
	double sum = 0;
	double sum_prog, sum2_prog;
	

//	Calcolo di r[i] per ogni blocco
				
	for(int i=0; i<N; i++){
		sum = 0;
		
		for(int j=0; j<L; j++)
			sum += rnd.Rannyu();
			
		r[i] = sum/L;
	}
	
//	Calcolo medie e incertezze progressive di r

	ofstream file("../output1_11.dat");


	for(int i=0; i<N; i++){
		sum_prog = 0;
		sum2_prog = 0;

		for(int j=0; j<=i; j++){
			sum_prog += r[j];
			sum2_prog += pow(r[j], 2);
		}
		sum_prog /= (i+1);
		sum2_prog /= (i+1);
		
		file << sum_prog-1./2 << " " << error(sum_prog, sum2_prog, i) << endl;
		
	} 
	
	file.close();	
	
	
///////////////////////	Parte 2 - Varianza ////////////////////	
	
//	Calcolo della varianza su ogni blocco
				
	for(int i=0; i<N; i++){
		sum = 0;
		
		for(int j=0; j<L; j++)
			sum += pow(rnd.Rannyu()-0.5, 2);
			
		sigma2[i] = sum/L;
	}
	
		
//	Calcolo medie e incertezze progressive della varianza

	file.open("../output1_12.dat");

	for(int i=0; i<N; i++){
		sum_prog = 0;
		sum2_prog = 0;

		for(int j=0; j<=i; j++){
			sum_prog += sigma2[j];
			sum2_prog += pow(sigma2[j], 2);
		}
		sum_prog /= (i+1);
		sum2_prog /= (i+1);
		
		file << sum_prog-1./12 << " " << error(sum_prog, sum2_prog, i) << endl;
		
	} 
	
	file.close();	
	
	
///////////////////////	Parte 3 - CHI Quadro //////////////////	
	
//	Osservazioni in ogni blocco
	double tmp;
	file.open("../output1_13.dat");	
				
	for(int i=0; i<N; i++){
		sum = 0;
		
		for(int j=0; j<L; j++){
			tmp = rnd.Rannyu();
			if(tmp>(i*1.)/N && tmp<(1.+i)/N) 
				sum++;
		}	
		file << pow(sum-100, 2)/(100) << endl;
	}

	file.close();

	


// Salvataggio del seme di uscita
   rnd.SaveSeed();
   return 0;
}

