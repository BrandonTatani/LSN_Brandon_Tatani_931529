/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include "Genetic_Algorithm.h"

using namespace std;

// -> Class City <-

	City :: City(int n, double r, Random* tmp){
		rnd = tmp;
		n_city = n;
		for(int i=0; i<n; i++){
			x.push_back(rnd->Rannyu(0., r));
			y.push_back(rnd->Rannyu(0., r));
		}
	};
	
	City :: ~City(){};
	
	double City :: GetL(int a, int b){
		return sqrt(pow(x[a]-x[b], 2) + pow(y[a]-y[b], 2));
	}
	
	int City :: GetN(){
		return n_city;
	}
	
	void City :: Print(){
		ofstream file("city.config");
		for(int i=0; i<n_city; i++)
			file << x[i] << " " << y[i] << endl;
		file.close();
	}
	
	double City :: GetX(int i){
		return x[i];
	}
	
	double City :: GetY(int i){
		return y[i];
	}
	
// -> Class Cromosoma <-

	Cromosoma :: Cromosoma(City *p, Random* tmp){
		rnd = tmp;
		citta = p;
		m_N = p-> GetN();
		m_L = 0;
		int gene;
		bool hit = false;
		vector <int> visita(m_N);
		for(int i=0; i<m_N; i++)
			visita[i] = 0;

		crom.push_back(1);
		for(int i=1; i<m_N; i++){
			hit = false;
			for(;hit==false;){
				gene = rnd->Rannyu(2, m_N);
				if(visita[gene-1] == 0){
					crom.push_back(gene);
					visita[gene-1] = 1;
					hit = true;
				}
			}
			
		}
		visita.clear();	
		UpdateL();	
	}
	
	Cromosoma :: Cromosoma(City* p, Random* tmp, Cromosoma a){
		rnd = tmp;
		citta = p;
		m_L = 0.;
		m_N = a.GetN();
		for(int i=0; i<m_N; i++)
			crom.push_back(a.Get_Gene(i));
	}
	
	Cromosoma :: ~Cromosoma(){}

	void Cromosoma :: Print(){
		cout << "[ ";
		for(int i=0; i<m_N; i++)
			cout << crom[i] << " ";
		cout << "]" << endl;
	}
	
	void Cromosoma :: UpdateL(){
		m_L = 0;
		for(int i=0; i<(m_N-1); i++)
			m_L += citta->GetL(crom[i]-1, crom[i+1]-1);
		m_L += citta -> GetL(crom[m_N-1]-1, 0);
		
	}
	
	double Cromosoma :: GetL(){
		return m_L;
	}
	
	int Cromosoma :: GetN(){
		return m_N;
	}
	
	void Cromosoma :: Set_Fit(double p){
		m_fit = p;
	}
	
	double Cromosoma :: Get_Fit(){
		return m_fit;
	}
	
	int Cromosoma :: Get_Gene(int a){
		return crom[a]; 
	}
	
	void Cromosoma :: Mutazione(double p){
		int m, n, a, b;
		if(rnd->Rannyu() < p){
			a = rnd->Rannyu(1,m_N-1);
			b = rnd->Rannyu(1,m_N-1);
			swap(crom[a], crom[b]);
		}	
		if(rnd->Rannyu() < p){
			vector <int> tmp(m_N);
			n = rnd->Rannyu(1,m_N-1);
			m = rnd->Rannyu(1,m_N-1);
			tmp[0] = 1;
			for(int i=0; i<m; i++)
				tmp[PBC(1+i+n,m_N)] = crom[PBC(1+i,m_N)];
			
			for(int i=0; i<m_N-m-1; i++)
				tmp[PBC(1+i+m+n,m_N)] = crom[PBC(1+m+i,m_N)];
			
			for(int i=1; i<m_N; i++)
				swap(crom[i], tmp[i]);
		}
		if(rnd->Rannyu() < p){
			m = rnd->Rannyu(1,int(m_N*0.5));
			for(int i=0; i<m; i++)
				swap(crom[1+i], crom[m_N-m+i]);				
		}
		if(rnd->Rannyu() < p){
			m = rnd->Rannyu(1,m_N-1);
			for(int i=0; i<m*0.5; i++)
				swap(crom[1+i], crom[m-i]);				
		}
		
	}
	
	void Cromosoma :: Crossover(Cromosoma a){
		int k=1;
		vector <int> visita(m_N);
		for(int i=0; i<m_N; i++){
			visita.push_back(0);
		}
		
		for(int i=0; i<m_N; i++){

			if(i<=floor(0.5*m_N))
				visita[crom[i]-1] = 1;
			else{
				while(visita[a.Get_Gene(k)-1] == 1){
					k++;
				}
					
				crom[i] = a.Get_Gene(k);
				visita[crom[i]-1] = 1;
				k++;
			}
		}	
			
	}
	
// -> Class Genetic <-

	Genetic :: Genetic(int n_cromosomi, City* citta, Random* rnd){
		popolo = n_cromosomi;
		m_rnd = rnd;
		m_citta = citta;
		p_mut = 0.2;
		p_cross = 0.7;
		m_p = 10;
		for(int i=0; i<popolo; i++){
			Cromosoma tmp(citta, rnd);
			membro.push_back(tmp);
		}
	}
	
	void Genetic :: Print(){
		for(int i=0; i<popolo; i++)
			membro[i].Print();
	}

	Genetic :: ~Genetic(){};
	
	void Genetic :: Print_Best(){
		m_file.open("Best.dat", ios::app);
//		m_file << membro[popolo-1].GetL() << endl;
		m_file << membro[0].GetL() << endl;
		m_file.close();
	}
	
	void Genetic :: Ordine(){
		Quicksort(membro, 0, membro.size()-1);
//		for(int i=0; i<popolo*0.5; i++)
//			swap(membro[i], membro[popolo-1-i]);
		
	}
	
	void Genetic :: Update_Param(){
		for(int i=0; i<popolo; i++){
//			cout << membro[i].GetL();
			membro[i].UpdateL();
		}
	}
	
	void Genetic :: Print_Param(){
//		cout << "Parametri dei cromosomi" << endl;
		for(int i=0; i<popolo; i++ )
			cout << membro[i].GetL() << " ";
		cout << endl;
	}
	
	void Genetic :: Mutazione(){
		for(int i=0; i<popolo; i++)
			membro[i].Mutazione(p_mut);
	}
	
	void Genetic :: Crossover(){
		vector <Cromosoma> figli;
		int gen_a, gen_b;
		for(int i=0; i<popolo*0.5; i++){
			gen_a = int(popolo * pow(m_rnd->Rannyu(), m_p));
			gen_b = int(popolo * pow(m_rnd->Rannyu(), m_p));
						
			Cromosoma a(membro[gen_a]);
			Cromosoma b(membro[gen_b]);

			if(m_rnd->Rannyu() < p_cross){
				a.Crossover(b);
				b.Crossover(a);
			}
			figli.push_back(a);
			figli.push_back(b);
		}
		for(int i=0; i<popolo; i++)
			swap(membro[i], figli[i]);
	}
	
	void Genetic :: Print_Cromosoma(){
		ofstream file("percorso.config");
		for(int i=0; i<m_citta->GetN(); i++){
			file << i << " " << m_citta->GetX(membro[0].Get_Gene(i)-1) << " " << m_citta->GetY(membro[0].Get_Gene(i)-1) << endl;
		}
		file.close();
	}
	
	void Genetic :: half_ave(){
		ofstream file("half_ave.dat", ios::app);
		double sum = 0;
		int n = popolo*0.5;
			for(int i=0; i<n; i++)
				sum += membro[i].GetL();
			file << sum/n << endl;
		file.close();
	}
// -> Funzioni Utili <- 


int Partition(vector <Cromosoma> &v, int start, int end){
	int pivot = end;
	int j = start;
	for(int i=start; i<end; i++){
		if(v[i].GetL() < v[pivot].GetL()){
			swap(v[i],v[j]);
			j++;
		}
	}
	swap(v[j], v[pivot]);
	return j;	
}

void Quicksort(vector <Cromosoma> &v, int start, int end){
	if(start < end){
		int p = Partition(v, start, end);
		Quicksort(v, start, p-1);
		Quicksort(v, p+1, end);
	}
}

int PBC(int i, int size){
	if(i >= size)
		return i-size+1;
	else
		return i;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
