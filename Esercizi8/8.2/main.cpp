#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include "random.h"

using namespace std;

struct Punti{
	double mu, sigma;
};
 
int main (int argc, char *argv[]){

// Creazione generatore di Random	

	Random rnd(0);

// Variabili dell'applicazione

	int M = 1e6, N=100, L=M/N;
	double av[N], sum_P, sum2_P;
	double a = 2.6;
	double x = 1.;
	ifstream input;
	ofstream output;
	ofstream config;
	
	double mu, sigma, tmp;
	
// Import parametri minimizzanti

	input.open("../8.1/config.min");
	input >> mu >> tmp >> sigma;
	input.close();
	

	
		
	rnd.SetV(mu, sigma);
	rnd.SetA(0);
	
	//	Calcolo di H[i] per ogni step
	
	config.open("config.final");

	for(int i=0; i<N; i++){
		av[i] = 0;
		for(int j=0; j<L; j++){
			rnd.Metropolis(x, a);
			config << x << endl;
			av[i] += rnd.Eval2(x);
		}
		av[i] /= double(L);
	}
	
	config.close();
		
	cout << "L'accettazione vale in media: " << rnd.GetA()/double(M) << endl;


//	Media a Blocchi e incertezza
	
	output.open("output8_3.dat");
	
	for(int i=1; i<=N; i++){
		sum_P = 0; 
		sum2_P = 0;
		for(int j=0; j<i; j++){

			sum_P += av[j]/i;
			sum2_P += pow(av[j],2)/i;			
		}

		output << i-1 << " " << sum_P << " " << error(sum_P, sum2_P, i-1) << endl;
	}	
	
	output.close();

// Salvataggio del seme di uscita
   rnd.SaveSeed();
   return 0;
}

