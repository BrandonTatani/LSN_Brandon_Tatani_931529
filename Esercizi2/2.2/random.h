

#ifndef __Random__
#define __Random__


class Random {

private:
  int m1,m2,m3,m4,l1,l2,l3,l4,n1,n2,n3,n4;

protected:

public:
  // constructors
  Random(int a);
  // destructor
  ~Random();
  // methods
  void SetRandom(int * , int, int);
  void SaveSeed();
  double Rannyu(void);
  double Rannyu(double min, double max);
  double Rannyu(double min, double max, double lambda);
  double Rannyu(double min, double max, double mean, double gamma);
  double Gauss(double mean, double sigma);
  double Exp(double lambda);
  double Reject();
  double Lorenziana(double mean, double gamma);
  double Angolo();
  int Walk();
  void Walk(double* x, double a);
};


double error(double AV, double AV2, int n);
double error2(double AV, double AV2, int n);
double error_rel(double AV, double AV2, int n);
double Coseno(double x);
double distanza(double* p);


#endif // __Random__


