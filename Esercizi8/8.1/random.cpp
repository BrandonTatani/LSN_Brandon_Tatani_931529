#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include "random.h"

using namespace std;


Random :: Random(int a){

	m_A=0;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
			input >> property;
			if( property == "RANDOMSEED" ){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				SetRandom(seed,p1,p2);
			}
      		}
      		input.close();
   	} else cerr << "PROBLEM: Unable to open seed.in" << endl;
}

Random :: ~Random(){}

void Random :: SaveSeed(){
   ofstream WriteSeed;
   WriteSeed.open("seed.out");
   if (WriteSeed.is_open()){
      WriteSeed << l1 << " " << l2 << " " << l3 << " " << l4 << endl;;
   } else cerr << "PROBLEM: Unable to open random.out" << endl;
  WriteSeed.close();
  return;
}

double Random :: Gauss(double mean, double sigma) {
   double s=Rannyu();
   double t=Rannyu();
   double x=sqrt(-2.*log(1.-s))*cos(2.*M_PI*t);
   return mean + x * sigma;
}

double Random :: Exp(double lambda){
	double p = Rannyu();
	return -log(1-p)/lambda;
}


double Random :: Lorenziana(double mean, double gamma){
	double p = Rannyu();
	return mean + gamma*tan(M_PI*(p-0.5));
}

double Random :: Angolo(){
	double x = Rannyu(-1, 1);
	double y = sqrt(1-x*x);
	return y;
}

double Random :: Reject(){
	double x, y;
	 
	do{
		x = Rannyu(0, 1);
		y = Rannyu(0,M_PI*0.5);
	}while(y >= 6./5.*(1-pow(x,2)*0.5));
	
	return x;
}

int Random :: Walk(){
	double passo = Rannyu();
	if(passo<1./3)
		return 0;
	else if(1./3<=passo && passo<2./3)
		return 1;
	else if(2./3<=passo)
		return 2;
}

void Random :: Walk(double* x, double a){
	double phi = Rannyu(0, M_PI);
	double theta = Rannyu(0, 2*M_PI);
		
	x[0] += a*sin(phi)*cos(theta);
	x[1] += a*sin(phi)*sin(theta);
	x[2] += a*cos(phi);
}

void Random :: Metropolis(double& x, double a){

	double p2 = x + Rannyu(-a, a);
		
	double A = min <double>(1, pow(Eval(p2), 2)/pow(Eval(x), 2));

	
	if(Rannyu() <= A){
		x = p2;
		m_A++;
	}
	
	return;
};
 
 //Rannyu tra Min e Max 
  
double Random :: Rannyu(double min, double max){
   return min+(max-min)*Rannyu();
}


double Random :: Rannyu(void){
  const double twom12=0.000244140625;
  int i1,i2,i3,i4;
  double r;

  i1 = l1*m4 + l2*m3 + l3*m2 + l4*m1 + n1;
  i2 = l2*m4 + l3*m3 + l4*m2 + n2;
  i3 = l3*m4 + l4*m3 + n3;
  i4 = l4*m4 + n4;
  l4 = i4%4096;
  i3 = i3 + i4/4096;
  l3 = i3%4096;
  i2 = i2 + i3/4096;
  l2 = i2%4096;
  l1 = (i1 + i2/4096)%4096;
  r=twom12*(l1+twom12*(l2+twom12*(l3+twom12*(l4))));

  return r;
}

void Random :: SetRandom(int * s, int p1, int p2){
  m1 = 502;
  m2 = 1521;
  m3 = 4071;
  m4 = 2107;
  l1 = s[0]%4096;
  l2 = s[1]%4096;
  l3 = s[2]%4096;
  l4 = s[3]%4096;
  l4 = 2*(l4/2)+1;
  n1 = 0;
  n2 = 0;
  n3 = p1;
  n4 = p2;

  return;
}

double Random :: GetA(){
	return m_A;
}

void Random :: SetA(int a){
	m_A = a;
}

double Random :: Eval(double x){
	return exp(-pow(x-mu,2)/(2.*sigma*sigma)) + exp(-pow(x+mu,2)/(2.*sigma*sigma));
}

double Random :: Eval2(double x){
	double piu = -pow((x-mu)/sigma, 2);
	double meno = -pow((x+mu)/sigma, 2);
	return (1./(2*sigma*sigma))*( (1+meno)*exp(0.5*meno) + (1+piu)*exp(0.5*piu))/Eval(x) + pot(x);
}

double pot(double x){
	return pow(x, 4) - 5.*0.5*pow(x, 2);
}


void Random :: SetV(double a, double b){
	mu = a;
	sigma = b;
}

double error(double AV, double AV2, int n){
	if(n==0)
		return 0;
	else
		return sqrt(fabs(AV2-pow(AV, 2))/n);
};

double error2(double AV, double AV2, int n){
	if(n==0)
		return 0;
	else if(AV2 == pow(AV, 2))
		return 0;
	else
		return sqrt(fabs(AV2-pow(AV, 2))/n)/(2*sqrt(AV));
};

double error_rel(double AV, double AV2, int n){
	if(n==0)
		return 0;
	else
		return sqrt((AV2-pow(AV, 2))/n)/AV;
};

double Coseno(double x){
	return M_PI*0.5 * cos(M_PI*0.5*x); 
}

double distanza(double* p){
	return sqrt(pow(p[0], 2) + pow(p[1], 2) + pow(p[2], 2));
}


void Print(double* p){
	cout << p[0] << " " << p[1] << " " << p[2] << endl;
	return;
}


/*double Max(double a, double b){
	if(a < b)
		return b;
	else
		return a;
};*/
