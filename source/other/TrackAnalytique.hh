#ifndef TrackAnalytique_h
#define TrackAnalytique_h 1

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

class TrackAnalytique
{
  public:
	TrackAnalytique();
	~TrackAnalytique();
	TrackAnalytique(vector<float> , vector<float> , vector<float> , vector<float> , vector<float> , vector<float> );
	TrackAnalytique(float* pointA1, float* pointA2);
	void Tracking();
	double* GetPointB1();
	double* GetPointB2();
	void Line(double t, double *p, double &x, double &y, double &z);
	TPolyLine3D* GetPolyLine();
	int GetSize();
	
	
	

  private:
	int size;
	//float* xa;
	//float* ya;
	//float* za;
	//float* errorXa;
	//float* errorYa;
	//float* errorZa;
  vector<float> Xa;//vecteur dans lequel sont stockes les positions x des hits
  vector<float> Ya;//vecteur dans lequel sont stockes les positions y des hits
  vector<float> Za;//vecteur dans lequel sont stockes les positions z des hits
  vector<float> errorXa;//error on the X position measurement, to define
  vector<float> errorYa;
  vector<float> errorZa;
  
	double* pointB1;
	double* pointB2;
	double parFit[4];

  
};



#endif