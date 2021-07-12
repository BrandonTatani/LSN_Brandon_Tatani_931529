#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <cstdlib>
using namespace std;

float error(vector <float> AV, vector <float> AV2, int n){
	if(n==0)
		return 0;
	else
		return sqrt((AV2[n]-pow(AV[n], 2))/n);
};

int main(){

	int M = 100000, N = 100, L = M/N;
	vector <float> rnd;
	ifstream input("random.dat");
	if(input.is_open()){

		float tmp = 0.;
		while(input >> tmp){
			rnd.push_back(tmp);
		}
		input.close();
		if(rnd.size() != M){
			cerr << "PROBLEM: random.dat does not contain the amount of numbers needed by the application:\t" << rnd.size() << 					" against " << M << " needed" << endl; 
		}
	}else{
		cerr << "PROBLEM: Unable to open the file containing randoms genereted" << endl;
	}

	vector <int> x;
	vector <float> av(N), av2(N), sum_prog(N), sum2_prog(N), err_prog(N);
	float sum = 0;
	
// Calcolo di r[i]
	
	for(int i=0; i<N; i++){
		sum = 0;
		for(int j=0; j<L; j++){
			sum += rnd[i*L+j];
		}
		av[i] = sum/L;
		av2[i] = pow(av[i], 2);
	}
	
// Calcolo delle incertezze

	for(int i=0; i<N; i++){
		for(int j=0; j<i+1; j++){
			sum_prog[i] += av[j];
			sum2_prog[i] += av2[j];
		}
		
		sum_prog[i]/=(i+1);
		sum2_prog[i]/=(i+1);
		err_prog[i] = error(sum_prog, sum2_prog, i);
	}	

// Calcolo dei lanci per blocco

	for(int i=0; i<N; i++)
		x.push_back(i*L);
		
// Stampa dei Risultati

	ofstream output("output.dat");
	for(int i=0; i<N; i++)
		output << x[i] << "\t" << sum_prog[i]-0.5 << "\t" << err_prog[i] << endl;
		
	output.close();
	

	return 0;
}
