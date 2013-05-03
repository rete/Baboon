#ifndef Track_h
#define Track_h 1

#include <TMath.h>
#include <TGraph2D.h>
#include <TRandom2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TH1.h>
#include <TVirtualFitter.h>
#include <TPolyLine3D.h>
#include <Math/Vector3D.h>
#include <vector>
#include <iostream>
using namespace std;

using namespace ROOT::Math;

class Track
{
  public:
	Track();
	~Track();
	//Track(int my_size, double* my_xa , double* my_ya, double* my_za);
	Track(vector<double> X, vector<double> Y, vector<double> Z);
	

	//void line(double t, double *p, double &x, double &y, double &z);
	//double distance2(double x,double y,double z, double *p);
	//void line();
	//double distance2();
	//void SumDistance2t(Int_t &, Double_t *, Double_t & sum, Double_t * par, Int_t ) ;
	void line(double t, double *p, double &x, double &y, double &z);
//	double distance2(double x,double y,double z, double *p);
	//void SumDistance2(int &, double *, double & sum, double * par, int );
	void Tracking();
	
	double* GetPointB1();
	double* GetPointB2();
	
	
	

  private:
	Int_t size;
	Double_t* xa;
	Double_t* ya;
	Double_t* za;
	//TGraph2D *gr;
	double t; 
	double* p; 
	double* x;
	double* y;
	double* z;
	double* pointB1;
	double* pointB2;
	//double* parFit;
	//bool first;
  
};


double distance2(double x,double y,double z, double *p);

void SumDistance2(int &, double *, double & sum, double * par, int );


#endif