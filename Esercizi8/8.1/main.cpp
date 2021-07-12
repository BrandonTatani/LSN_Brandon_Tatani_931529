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

// Controllo sui parametri

	if(argc != 5)
		cerr << "Uso del programma: ./main sigma_min sigma_max mu_min mu_max" << endl;

// Creazione generatore di Random	

	Random rnd(0);

// Variabili dell'applicazione

	int M = 1e6, N=100, L=M/N;
	double av[N], sum_P, sum2_P;
	double a = 2.6;
	double x = 1.;
	const int set = 100;
	Punti V[set];
	rnd.SetV(1.,1.);
	ofstream file;
	
	double min = 0.;
	int tmp = 0;
	
//	Set della griglia (mu, sigma)

	double mu_min = atof(argv[3]), mu_max = atof(argv[4]);
	double  sigma_min = atof(argv[1]), sigma_max = atof(argv[2]);
	
	double d_mu = (mu_max-mu_min)/10.;
	double d_sigma = (sigma_max-sigma_min)/10.;

	for(int i=0; i<10; i++){
		for(int k=0; k<10; k++){
			V[i*10+k].mu = mu_min + i*d_mu;
			V[i*10+k].sigma = sigma_min + k*d_sigma;
		}
	}
	
	file.open("output8_2.dat");
	
	for(int t=0; t<set; t++){
		
		rnd.SetV(V[t].mu, V[t].sigma);
		rnd.SetA(0);
		
		//	Calcolo di H[i] per ogni step

		for(int i=0; i<N; i++){
			av[i] = 0;
			for(int j=0; j<L; j++){
				rnd.Metropolis(x, a);
				av[i] += rnd.Eval2(x);
			}
			av[i] /= double(L);
		}
		

	//	Media a Blocchi e incertezza
		
		for(int i=1; i<=N; i++){
			sum_P = 0; 
			sum2_P = 0;
			for(int j=0; j<i; j++){
	
				sum_P += av[j]/i;
				sum2_P += pow(av[j],2)/i;			
			}
			if(i==N){
				file << V[t].mu << " " << V[t].sigma << " " << sum_P << " " << error(sum_P, sum2_P, i-1) << endl;
				if(t==0){
					min = sum_P;
					tmp = t;
				}
				if(sum_P < min){
					min = sum_P;
					tmp = t;
				}
			}
		}
		
		if(t%10 == 0){
			cout << t << endl;
			cout << "L'accettazione vale in media: " << rnd.GetA()/double(M) << endl;
		}
	}
	
	file.close();
	
	//	Print parametri che minimizzano H
	
	file.open("config.min");
	
	file << V[tmp].mu << " " << d_mu << " " << V[tmp].sigma << " " << d_sigma << endl;
	
	file.close();

// Salvataggio del seme di uscita
   rnd.SaveSeed();
   return 0;
}

