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
		City(int n);
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
		~Cromosoma();

		void UpdateL();	
		double GetL();
		void Set_Fit(double);
		double Get_Fit();
		void Print();
		void Mutazione(double);
		void Crossover(Cromosoma);
		int Get_Gene(int);
		void Set_Gene(int, int); 
		int GetN();
		int* Get_Address(int);
		
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
		void Print(int);		
		void Print_Param();
		int Get_Gene(int, int); 
		void Set_Gene(int, int, int);
		int* Get_Address(int, int);
		void Copy(int*, int);
		double GetTmp(){return m_tmp;};
		double Get_PCross();
		void Set_PCross(double);
		void Print_Best();
		void Print_Cromosoma();
		void half_ave();
		double GetL(int);		
				
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

int Partition(vector <Cromosoma> &v, int start, int end);
void Quicksort(vector <Cromosoma> &v, int start, int end);
int PBC(int, int);

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
