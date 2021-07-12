#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

// Creazione generatore di Random	
  
	int seed[4] = {0, 0, 1, 1} ;
	Random rnd(seed);

// Variabili dell'applicazione

	int M = 10000, N=100, n=M/N;
	double *p = new double(3);
	double d[N], d2[N], tmp[N], tmp2[N];
	
	double sum_prog, sum2_prog;

	for(int i=0; i<N; i++){
		d[i]=0;
		d2[i]=0;
	}

//	Calcolo di f[i] per ogni blocco

	for(int i=0; i<N; i++)
		tmp[i] = 0;
	
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
		}

	}
	
	ofstream file("../output2_21.csv");
	
	file << "RMS,Errore" << endl;
	
	for(int i=0; i<N; i++){
		file << sqrt(d[i]) << "," << error2(d[i], d2[i], N) << endl;
	}
	
	file.close();
				
/*	for(int i=0; i<M; i++){
		p[0] = 0;
		p[1] = 0;
		p[2] = 0;
		
		for(int j=1; j<N; j++){
			if(rnd.Rannyu() < 0.5)
				p[rnd.Walk()]--;
			else
				p[rnd.Walk()]++;
			
			d[j] += pow(distanza(p), 2)/M;
		}

	}

//	Calcolo medie e incertezze progressive di I

	ofstream file("../output2_21.csv");

	file << "Distanza" << "," << "errore" << endl;

	for(int i=0; i<N; i++){
		sum_prog = 0;
		sum2_prog = 0;

		for(int j=0; j<=i; j++){
			sum_prog += d[j];
			sum2_prog += pow(d[j], 2);
		}
		sum_prog /= (i+1);
		sum2_prog /= (i+1);
		
	//	file << sqrt(sum_prog) << "," << error2(sum_prog, sum2_prog, i) << endl;
		file << sqrt(d[i]) << "," << error2(sum_prog, sum2_prog, i) << endl;
	} 
	
	file.close();	
*/	

// 							-> Parte Due <- 						//

	for(int i=0; i<N; i++){
		d[i]=0;
		d2[i]=0;
	}

//	Calcolo di f[i] per ogni blocco
				
	for(int i=0; i<M; i++){
		p[0] = 0;
		p[1] = 0;
		p[2] = 0;
		
		for(int j=1; j<N; j++){
			rnd.Walk(p, 1.);
			d[j] += pow(distanza(p), 2)/M;
			d2[j] += pow(distanza(p), 4)/M;
		}

	}

//	Calcolo medie e incertezze progressive di I

	file.open("../output2_22.csv");

	file << "Distanza" << "," << "errore" << endl;

	for(int i=0; i<N; i++){
		sum_prog = 0;
		sum2_prog = 0;

		for(int j=0; j<=i; j++){
			sum_prog += d[j];
			sum2_prog += d2[j];
		}
		sum_prog /= (i+1);
		sum2_prog /= (i+1);
		
	//	file << sqrt(sum_prog) << "," << error2(sum_prog, sum2_prog, i) << endl;
	//	file << sqrt(sum_prog) << "," << error(sum_prog, sum2_prog, i)*sqrt(sum_prog) << endl;
		file << sqrt(d[i]) << "," << error2(d[i], d2[i], i) << endl;
	} 
	
	file.close();	
	

	
// Salvataggio del seme di uscita
   rnd.SaveSeed();
   return 0;
}

