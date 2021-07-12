#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

// Creazione generatore di Random	
  
	Random rnd(a);

// Variabili dell'applicazione

	double sum = 0;
	int M = 10000;
	int lambda = 1;
	int mean=0, gamma=1;
	int N[4] = {1, 2, 10, 100};
	
// Figura 1 -> Distribuzione uniforme
	
	ofstream file("../figura_1.csv");
	
	file << "N1,N2,N10,N100" << endl;
	
	for(int i=0; i<M; i++){
		file << i << ",";
		for(int j=0; j<4; j++){
			sum=0;
			for(int k=0; k<N[j]; k++){
				sum += rnd.Rannyu()/N[j];
			}
			if(j==3)	file << sum << endl;
			else	file << sum << ",";
		}
		
	}
	
	file.close();


// Figura 2 -> Distribuzione esponenziale

	file.open("../figura_2.csv");
	
	file << "N1,N2,N10,N100" << endl;
	
	for(int i=0; i<M; i++){
		file << i << ",";
		for(int j=0; j<4; j++){
			sum=0;
			for(int k=0; k<N[j]; k++){
				sum += rnd.Exp(lambda)/N[j];
			}
			if(j==3)	file << sum << endl;
			else	file << sum << ",";
		}
		
	}
	
	file.close();

// Figura 3 -> Distribuzione Lorenziana
	
	file.open("../figura_3.csv");
	
	file << "N1,N2,N10,N100" << endl;
	
	for(int i=0; i<M; i++){
		file << i << ",";
		for(int j=0; j<4; j++){
			sum=0;
			for(int k=0; k<N[j]; k++){
				sum += rnd.Lorenziana(mean, gamma)/N[j];
			}
			if(j==3)	file << sum << endl;
			else	file << sum << ",";
		}
		
	}
	
	file.close();

// Salvataggio del seme di uscita
   rnd.SaveSeed();
   return 0;
}

