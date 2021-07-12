/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#ifndef __GENE__
#define __GENE__

#include <iostream>
#include <fstream>
#include "random.h"
#include <vector>
#include <algorithm>
using namespace std;

//Città:

class City {

	public:
		
		City(int, double, Random*);
		~City();
	
		double GetL(int, int); // L^2 tra la città A e la città B
		int GetN(); // Recuperà il numero di città
		void Print();
		double GetX(int);
		double GetY(int);
	
	private:
	
		vector <double> x;
		vector <double> y;
		int n_city;
		Random* rnd;

};

class Cromosoma {

	public:
		
		Cromosoma(City*, Random*);
		Cromosoma(City*, Random*, Cromosoma);
		Cromosoma(City*, Random*, Cromosoma*);
		~Cromosoma();

		void UpdateL();	
		double GetL();
		void Set_Fit(double);
		double Get_Fit();
		void Print();
		void Mutazione(double);
		void Mutazione(int);
		void Crossover(Cromosoma);
		
		void Swap(Cromosoma*);
		void Swap(Cromosoma);
		int Get_Gene(int);
		int GetN();
		
	private:
		
		City* citta;
		Random* rnd;
		int m_N;
		double m_L; // Lunghezza totale del tratto;
		double m_fit; // Probabilità di essere scelto come Genitore;
		vector <int> crom;
	
};

class Genetic {

	public: 
	
		Genetic(int, City*, Random*);
		~Genetic();

		void Ordine();
		void Update_Param();
		void Crossover();
		void Mutazione();
		
		void Print();		
		void Print_Param();
		
		double GetTmp(){return m_tmp;};
		double Get_PCross();
		void Set_PCross(double);
		void Print_Best();
		void Print_Cromosoma();
				
	private:
	
		int popolo;
		ofstream m_file;
		vector <Cromosoma> membro;
//		vector <double> param;
//		vector <double> fitness;
		Random* m_rnd;
		City* m_citta;
		double p_mut, p_cross, m_p, m_tmp=0;
};

class Annealing {
	public:
		Annealing(City*, Random*, double);
		~Annealing();
		
		
		void Update_Param();
		void Step();
		void Print();		
		void Print_Param();
		void Save_Param();
		void Save_Percorso();
		double Boltzmann(double);


		double GetT();
		void SetT(double);

	private:
		Cromosoma* membro;
		Random* m_rnd;
		City* m_citta;
		ofstream m_file;
		double p_mut, m_T, beta;
};
int Partition(vector <Cromosoma> &v, int start, int end);
void Quicksort(vector <Cromosoma> &v, int start, int end);
int PBC(int, int);
int i_max(vector <int>);
#endif

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
