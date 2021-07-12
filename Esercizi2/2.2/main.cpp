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

	int M = 10000, N=100, n=M/N;
	double *p = new double(3);
	double d[N], d2[N], tmp[N], tmp2[N];


	for(int i=0; i<N; i++){
		d[i]=0;
		d2[i]=0;
		tmp[i] = 0;
		tmp2[i] = 0;
	}

//	Calcolo di r[i] per ogni step


	for(int i=0; i<N; i++){
		for(int j=0; j<n; j++){
			p[0] = 0;
			p[1] = 0;
			p[2] = 0;
			
			for(int k=1; k<N; k++){
				if(rnd.Rannyu() < 0.5)
					p[rnd.Walk()]--;
				else
					p[rnd.Walk()]++;
				
				tmp[k] += pow(distanza(p), 2)/n;
				tmp2[k] += pow(distanza(p), 4)/n;
			}
		}

		for(int s=0; s<N; s++){
			d[s] += tmp[s]/N;
			d2[s] += tmp2[s]/N;
			tmp[s] = 0;
			tmp2[s] = 0;
		}

	}
	
//	Stampa dei risultati
	
	ofstream file("../output2_21.csv");
	
	file << "RMS,Errore" << endl;
	
	for(int i=0; i<N; i++){
		file << sqrt(d[i]) << "," << error2(d[i], d2[i], N) << endl;
	}
	
	file.close();
				


// 							-> Parte Due <- 						//

	for(int i=0; i<N; i++){
		d[i]=0;
		d2[i]=0;
		tmp[i] = 0;
		tmp2[i] = 0;
	}

//	Calcolo di r[i] per ogni step

	
	for(int i=0; i<N; i++){
		for(int j=0; j<n; j++){
			p[0] = 0;
			p[1] = 0;
			p[2] = 0;
			
			for(int k=1; k<N; k++){
				rnd.Walk(p, 1);
				
				tmp[k] += pow(distanza(p), 2)/n;
				tmp2[k] += pow(distanza(p), 4)/n;
			}
		}

		for(int s=0; s<N; s++){
			d[s] += tmp[s]/N;
			d2[s] += tmp2[s]/N;
			tmp[s] = 0;
			tmp2[s] = 0;
		}

	}
	
//	Stampa dei risultati
	
	file.open("../output2_22.csv");
	
	file << "RMS,Errore" << endl;
	
	for(int i=0; i<N; i++){
		file << sqrt(d[i]) << "," << error2(d[i], d2[i], N) << endl;
	}
	
	file.close();

	
// Salvataggio del seme di uscita
   rnd.SaveSeed();
   return 0;
}

