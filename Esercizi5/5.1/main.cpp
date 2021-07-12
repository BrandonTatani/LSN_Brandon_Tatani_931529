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

	int M = 1e6, N=100, L=M/N;
	double av[N], sum_P, sum2_P;
	double a = 0.7;
	double* p1 = new double(3);
	double* p2 = new double(3);

	
	for(int i=0; i<3; i++){
		p1[i] = 1.;
		p2[i] = 1.;
	}
	
//	Calcolo di r[i] per ogni step

	for(int i=0; i<N; i++){
		av[i] = 0;
		for(int j=0; j<L; j++){
			rnd.Metropolis(p1, p2, a, 1);
			av[i] += distanza(p1)/double(L);
		}
	}
	
	cout << "L'accettazione vale in media: " << rnd.GetA()/double(M) << endl;


//	Media a Blocchi e incertezza

	ofstream file("../output5_5.dat");
	
	for(int i=1; i<=N; i++){
		sum_P = 0; 
		sum2_P = 0;
		for(int j=0; j<i; j++){

			sum_P += av[j]/i;
			sum2_P += pow(av[j],2)/i;			
		}
		file << sum_P << " " << error(sum_P, sum2_P, i-1) << endl;
	}

	file.close();

	
//		-> Funzione d'onda 2,1,0 <-

	for(int i=0; i<3; i++){
		p1[i] = 1.;
		p2[i] = 1.;
	}

	rnd.SetA(0);

//	Calcolo di r[i] per ogni step
	a = 1.8;
	
	for(int i=0; i<N; i++){
		av[i] = 0;
		for(int j=0; j<L; j++){
			rnd.Metropolis(p1, p2, a, 2);
			av[i] += distanza(p1)/double(L);
		}
	}
	
	cout << "L'accettazione vale in media: " << rnd.GetA()/double(M) << endl;


//	Media a Blocchi e incertezza

	file.open("../output5_6.dat");
	
	for(int i=1; i<=N; i++){
		sum_P = 0; 
		sum2_P = 0;
		for(int j=0; j<i; j++){

			sum_P += av[j]/i;
			sum2_P += pow(av[j],2)/i;			
		}
		file << sum_P << " " << error(sum_P, sum2_P, i-1) << endl;
	}

	file.close();


/*
//		-> Simulazione partendo da un punto molto lontano <-
	
// Punto di partenza (x,y,z) = (10,10,10)

	for(int i=0; i<3; i++){
		p1[i] = 10.;
		p2[i] = 10.;
	}

//	Calcolo di r[i] per ogni step

	M = 1e3;
	L = M/N;
	a = 2.5;

	for(int i=0; i<N; i++){
		av[i] = 0;
		for(int j=0; j<L; j++){
			rnd.Metropolis(p1, p2, a, 1);
			av[i] += distanza(p1)/double(L);
		}
	}
	
//	Media a Blocchi e incertezza

	file.open("../output5_3.dat");
	
	for(int i=1; i<=N; i++){
		sum_P = 0; 
		sum2_P = 0;
		for(int j=0; j<i; j++){

			sum_P += av[j]/i;
			sum2_P += pow(av[j],2)/i;			
		}
		file << sum_P << " " << error(sum_P, sum2_P, i-1) << endl;
	}

	file.close();
	
// Punto di partenza (x,y,z) = (1,1,1)

	for(int i=0; i<3; i++){
		p1[i] = 1.;
		p2[i] = 1.;
	}

//	Calcolo di r[i] per ogni step

	for(int i=0; i<N; i++){
		av[i] = 0;
		for(int j=0; j<L; j++){
			rnd.Metropolis(p1, p2, a, 1);
			av[i] += distanza(p1)/double(L);
		}
	}
	
//	Media a Blocchi e incertezza

	file.open("../output5_4.dat");
	
	for(int i=1; i<=N; i++){
		sum_P = 0; 
		sum2_P = 0;
		for(int j=0; j<i; j++){

			sum_P += av[j]/i;
			sum2_P += pow(av[j],2)/i;			
		}
		file << sum_P << " " << error(sum_P, sum2_P, i-1) << endl;
	}

	file.close();
*/


// Salvataggio del seme di uscita
   rnd.SaveSeed();
   return 0;
}

