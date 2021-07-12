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
	double C[N], C2[N], P[N], P2[N],tmp, tmp2, sum_P, sum2_P;
	double S_T = 0, S_0=100, K=100, T=1., r=0.1, sigma=0.25;


	for(int i=0; i<N; i++){
		C[i]=0;
		C2[i]=0;
		P[i]=0;
		P2[i]=0;
	}

	
//	Calcolo di C[i] per ogni step

	for(int i=0; i<N; i++){
		for(int j=0; j<n; j++){
			S_T = S_0*exp((r-0.5*pow(sigma, 2))*T + sigma*rnd.Gauss(0,T));
			C[i] += exp(-r*T)*max <double>(0, S_T-K)/n;
			P[i] += exp(-r*T)*max <double>(0, K-S_T)/n;
	
		}
		C2[i] = pow(C[i], 2);
		P2[i] = pow(P[i], 2);
	}

//	Media a Blocchi e incertezza

	ofstream file("../output3_1.csv");
	
	file << "Call,Err_Call,Put,Err_Put" << endl;

	for(int i=1; i<=N; i++){
		tmp = 0;
		tmp2 = 0;
		sum_P = 0; 
		sum2_P = 0;
		for(int j=0; j<i; j++){
			tmp += C[j]/i;
			tmp2 += C2[j]/i;
			sum_P += P[j]/i;
			sum2_P += P2[j]/i;			
		}
		file << tmp << "," << error(tmp, tmp2, i-1) << "," << sum_P << "," << error(sum_P, sum2_P, i-1) << endl;
	}

	file.close();
	
	
//		-> Call Price with S(t) simulation <-


//	Calcolo di C[i] per ogni step

	for(int i=0; i<N; i++){
		C[i] = 0;
		C2[i] = 0;
		P[i] = 0;
		P2[i] = 0;
		for(int j=0; j<n; j++){
			S_0 = 100.;
			for(double t=0; t<T; t+=T/N)
				S_0 *= exp((r-0.5*pow(sigma,2))*(T/100.) + sigma*rnd.Gauss(0,1)*sqrt(T/100.));
			
			C[i] += exp(-r*T)*max <double>(0, S_0-K)/n;
			P[i] += exp(-r*T)*max <double>(0, K-S_0)/n;
		}
		C2[i] = pow(C[i], 2);
		P2[i] = pow(P[i], 2);
	}

//	Media a Blocchi e incertezza

	file.open("../output3_2.csv");
	
	file << "Call,Err_Call,Put,Err_Put" << endl;

	for(int i=1; i<=N; i++){
		tmp = 0;
		tmp2 = 0;
		sum_P = 0;
		sum2_P = 0;
		for(int j=0; j<i; j++){
			tmp += C[j]/i;
			tmp2 += C2[j]/i;
			sum_P += P[j]/i;
			sum2_P += P2[j]/i;			
		}
		file << tmp << "," << error(tmp, tmp2, i-1) << "," << sum_P << "," << error(sum_P, sum2_P, i-1)<< endl;
	}

	file.close();
	

// Salvataggio del seme di uscita
   rnd.SaveSeed();
   return 0;
}

