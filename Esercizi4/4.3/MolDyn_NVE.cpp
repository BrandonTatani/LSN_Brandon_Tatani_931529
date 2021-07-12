/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(int argc, char** argv){ 
	if(argc!=2 || (atoi(argv[1])!=0 && atoi(argv[1])!=1)){
		cerr << "Use of the program: <./main.exe 0> or <./main.exe 1> based on the existance of the file 'old.0'\n i.e. (0) if old.0 doesn't exists\t (1) if old.0 exists" << endl;
		return 2;
		}
	bool old_exist = atoi(argv[1]);
  Input(old_exist);             //Inizialization
  int nconf = 1;
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();     //Properties measurement
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 

		
        nconf += 1;
     }	
     if(istep == nstep-1)
	ConfOld();
  }
  ConfFinal();         //Write final configuration to restart

	Blocking(nstep/10);

  return 0;
}

/***********************************************************************************/

void Input(bool old_exist){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();
  
  if(old_exist == false){

	//Prepare initial velocities
 	  cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   	double sumv[3] = {0.0, 0.0, 0.0};
   	for (int i=0; i<npart; ++i){
   	  vx[i] = rand()/double(RAND_MAX) - 0.5;
   	  vy[i] = rand()/double(RAND_MAX) - 0.5;
   	  vz[i] = rand()/double(RAND_MAX) - 0.5;
	
   	  sumv[0] += vx[i];
   	  sumv[1] += vy[i];
   	  sumv[2] += vz[i];
   	}
   	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   	double sumv2 = 0.0, fs;
   	for (int i=0; i<npart; ++i){
   	  vx[i] = vx[i] - sumv[0];
   	  vy[i] = vy[i] - sumv[1];
   	  vz[i] = vz[i] - sumv[2];

   	  sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   	}
   	sumv2 /= (double)npart;

   	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	for (int i=0; i<npart; ++i){
		vx[i] *= fs;
  	   	vy[i] *= fs;
     		vz[i] *= fs;
	
   		xold[i] = Pbc(x[i] - vx[i] * delta);
    		yold[i] = Pbc(y[i] - vy[i] * delta);
    		zold[i] = Pbc(z[i] - vz[i] * delta);
  	}
 		return;
	}else{

  cout << "Read configuration in (t-dt) from file old.0 " << endl << endl;
  ReadConf.open("old.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> xold[i] >> yold[i] >> zold[i];
    xold[i] = xold[i] * box;
    yold[i] = yold[i] * box;
    zold[i] = zold[i] * box;
  }
  ReadConf.close();
  
// First step

	Move();

// Temperature estimation
	double stima_temp = 0;
	double t = 0;

	 for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    cout << "temp = " << stima_temp << endl;
    
// Fattore di scala
	double fs = sqrt(temp/stima_temp);

// Riscalamento
/*	for (int i=0; i<npart; ++i){
	     vx[i] *= fs;
	     vy[i] *= fs;
	     vz[i] *= fs;

	     xold[i] = Pbc(x[i] - vx[i] * delta);
	     yold[i] = Pbc(y[i] - vy[i] * delta);
	     zold[i] = Pbc(z[i] - vz[i] * delta);
	}	
*/
   return;
}

return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("output_epot.dat",ios::app);
  Ekin.open("output_ekin.dat",ios::app);
  Temp.open("output_temp.dat",ios::app);
  Etot.open("output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);

//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; //Potential energy per particle
    stima_kin = t/(double)npart; //Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; //Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfOld(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print semi-final configuration to file config_semi.final " << endl << endl;
  WriteConf.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double error(double av, double av2, int i){
	if(i==0){
		return 0;
	}else{
		return sqrt(fabs(av-av2)/(i));
	}
}

void Blocking(int nstep){
	const int N = 100;
	int L = nstep/N;
	ifstream f_Temp, f_Etot, f_Ekin, f_Epot;
	double av_T[N], av_E[N], av_K[N], av_V[N];
	double sum_prog[4], sum_prog2[4], tmp=0;
	
	/*
		0 -> temp
		1 -> Etot
		2 -> Ekin
		3 -> Epot
	*/
	
	f_Temp.open("output_temp.dat", ios::in);
	f_Etot.open("output_etot.dat", ios::in);
	f_Ekin.open("output_ekin.dat", ios::in);
	f_Epot.open("output_epot.dat", ios::in);
		
	for(int i=0; i<N; i++){
		
		av_T[i] = 0;
		av_E[i] = 0;
		av_K[i] = 0;
		av_V[i] = 0;	
			
		for(int j=0; j<L; j++){
			f_Temp >> tmp;
			av_T[i] += tmp/L;
			f_Etot >> tmp;
			av_E[i] += tmp/L;
			f_Ekin >> tmp;
			av_K[i] += tmp/L;
			f_Epot >> tmp;
			av_V[i] += tmp/L;

		}
	}
	
	f_Temp.close();
	f_Etot.close();
	f_Ekin.close();
	f_Epot.close();
	
	// Blocchi
	
	ofstream ave_temp, ave_etot, ave_ekin, ave_epot;
	ave_temp.open("ave_temp.dat");
	ave_etot.open("ave_etot.dat");
	ave_ekin.open("ave_ekin.dat");
	ave_epot.open("ave_epot.dat");
	
	
	for(int i=0; i<N; i++){
		for(int j=0; j<4; j++){
			sum_prog[j]=0;
			sum_prog2[j]=0;
		}
		for(int j=0; j<=i; j++){
			sum_prog[0] += av_T[j]/(i+1);
			sum_prog2[0] += pow(av_T[j], 2)/(i+1);
			sum_prog[1] += av_E[j]/(i+1);
			sum_prog2[1] += pow(av_E[j], 2)/(i+1);
			sum_prog[2] += av_K[j]/(i+1);
			sum_prog2[2] += pow(av_K[j], 2)/(i+1);
			sum_prog[3] += av_V[j]/(i+1);
			sum_prog2[3] += pow(av_V[j], 2)/(i+1);
		}
		
		ave_temp << sum_prog[0] << " " << error(sum_prog[0], sum_prog2[0], i) << endl;
		ave_etot << sum_prog[1] << " " << error(sum_prog[1], sum_prog2[1], i) << endl;
		ave_ekin << sum_prog[2] << " " << error(sum_prog[2], sum_prog2[2], i) << endl;
		ave_epot << sum_prog[3] << " " << error(sum_prog[3], sum_prog2[3], i) << endl;
	}
	
	ave_temp.close();
	ave_etot.close();
	ave_ekin.close();
	ave_epot.close();
	
	return;
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
