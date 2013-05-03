#include "TrackAnalytique.hh"



TrackAnalytique::TrackAnalytique()
{
  //cout<<"----- in track analytique"<<endl;
 	pointB1 = new double[3];
	pointB2 = new double[3]; 
	pointB1[0] = 0;
	pointB1[1] = 0;
	pointB1[2] = 0;
	pointB2[0] = 1;
	pointB2[1] = 1;
	pointB2[2] = 1;	
}

TrackAnalytique::TrackAnalytique(vector<float> X, vector<float> Y, vector<float> Z, vector<float> errorX, vector<float> errorY, vector<float> errorZ)
//Track::Track(vector<Double_t> X, vector<Double_t> Y, vector<Double_t> Z)
{
  //cout<<"---- here 1 ----"<<endl;
  	pointB1 = new double[3];
	//cout<<"---- here 2 ----"<<endl;
	pointB2 = new double[3]; 
	/*pointB1[0] = -1;
	pointB1[1] = 1;
	pointB1[2] = 3;
	pointB2[0] = -2;
	pointB2[1] = 2;
	pointB2[2] = 4;*/
	//pointB1[0] = 0;
	//pointB1[1] = 0;
	//pointB1[2] = 0;
	//pointB2[0] = 1;
	//pointB2[1] = 1;
	//pointB2[2] = 1;	

  //cout<<"---- here 3 ----"<<endl;
size = X.size();
Xa = X;
Ya = Y;
Za = Z;
errorXa = errorX;
errorYa = errorY;
errorZa = errorZ;
  
  
  //cout<<"----------Construction de l'objet track OK----------"<<endl;
}


TrackAnalytique::TrackAnalytique(float* pointA1, float* pointA2)
{
  //cout<<"---- constructor first track"<<endl;
  	pointB1 = new double[3];
	pointB2 = new double[3]; 
	//pointB1[0] = -1;
	//pointB1[1] = 1;
	//pointB1[2] = 3;
	//pointB2[0] = -2;
	//pointB2[1] = 2;
	//pointB2[2] = 4;	
  
  Xa.push_back(pointA1[0]);
  Xa.push_back(pointA2[0]);
  Ya.push_back(pointA1[1]);
  Ya.push_back(pointA2[1]);
  Za.push_back(pointA1[2]);
  Za.push_back(pointA2[2]);
  
  errorXa.push_back(1);
  errorXa.push_back(1);
  errorYa.push_back(1);
  errorYa.push_back(1);
  errorZa.push_back(1);
  errorZa.push_back(1);
  
  //cout<<"----------Construction de l'objet track OK----------"<<endl;
}

TrackAnalytique::~TrackAnalytique()
{
	//delete [] xa;
	//delete [] ya;
	//delete [] za;
	//delete [] errorXa;
	//delete [] errorYa;
	//delete [] errorZa;
	//delete [] pointB1;
	//delete [] pointB2;
		  Xa.clear();
		  Ya.clear();
		  Za.clear();
		  errorXa.clear();
		  errorYa.clear();
		  errorZa.clear();
  //cout<<"----------Passage destructeur track OK----------"<<endl;
}

/*
void TrackAnalytique::Tracking()
{
 
    if(size==2)
	{
		pointB1[0] = xa[0];
		pointB1[1] = ya[0];
		pointB1[2] = za[0];
		pointB2[0] = xa[1];
		pointB2[1] = ya[1];
		pointB2[2] = za[1];
	}
	
	if(size>2)
	{
  //-------------------oooooOOOOO00000OOOOOooooo---------------------#
  //                                            
  //	THE PARAMETRIC LINE EQUATION	  
  //                                       
  //-------------------oooooOOOOO00000OOOOOooooo---------------------#

  // a parameteric line is define from 6 parameters but 4 are independent
  // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
  // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1; 
  // x = p[0] + p[1]*t; 
  // y = p[2] + p[3]*t;
  // z = t;  


  //-------------------oooooOOOOO00000OOOOOooooo---------------------#
  //                                            
  //	EQUATIONS TO SOLVE WITH "moindres carres" METHOD
  //                                       
  //-------------------oooooOOOOO00000OOOOOooooo---------------------#

  // x = p[0] + p[1]*z; -> equation 1
  // y = p[2] + p[3]*z; -> equation 2
  
//Moindres Carres ordinaires : non prise en compte du poids 
double xsum, ysum, zsum, zzsum, xzsum, yzsum;

xsum = 0.0;
ysum = 0.0;
zsum = 0.0;
zzsum = 0.0;
xzsum = 0.0;
yzsum = 0.0;


for (int i=0;i<size;i++)
{
  //for equation 1
  zsum = zsum + za[i];
  xsum = xsum + xa[i];   
  zzsum = zzsum + (za[i]*za[i]);
  xzsum = xzsum + xa[i]*za[i];
  
  
  //for equation 2
  ysum = ysum + ya[i];  
  yzsum = yzsum + ya[i]*za[i]/(errorYa[i]*errorYa[i]);  
  
cout<<setw(15)<<"zsum"<<setw(15)<<"xsum"<<setw(15)<<"zzsum"<<setw(15)<<"xzsum"<<setw(15)<<"ysum"<<setw(15)<<"yzsum"<<endl;
cout<<setw(15)<<zsum<<setw(15)<<xsum<<setw(15)<<zzsum<<setw(15)<<xzsum<<setw(15)<<ysum<<setw(15)<<yzsum<<endl<<endl;
  
}

double A1 = zsum;
double B1 = size;
double C1 = xsum;
double D1 = zzsum;
double E1 = xzsum;

double C2 = ysum;
double E2 = yzsum;


cout<<setw(15)<<"A1"<<setw(15)<<"B1"<<setw(15)<<"C1"<<setw(15)<<"D1"<<setw(15)<<"E1"<<setw(15)<<"C2"<<setw(15)<<"E2"<<endl;
cout<<setw(15)<<A1<<setw(15)<<B1<<setw(15)<<C1<<setw(15)<<D1<<setw(15)<<E1<<setw(15)<<C2<<setw(15)<<E2<<endl<<endl;


double parFit[4];
parFit[0] = (D1*C1-E1*A1)/(B1*D1-A1*A1);
parFit[1] = (E1*B1-C1*A1)/(B1*D1-A1*A1);
parFit[2] = (D1*C2-E2*A1)/(B1*D1-A1*A1);
parFit[3] = (E2*B1-C2*A1)/(B1*D1-A1*A1);

for(int i=0; i<4; i++)
{
  cout<<"parFit["<<i<<"] : "<<parFit[i]<<endl;
}



  //-------------------oooooOOOOO00000OOOOOooooo---------------------#
  //                                            
  //	FITTED LINE'S INTERSECTION WITH THE FIRST AND THE LAST PLANES
  //                                       
  //-------------------oooooOOOOO00000OOOOOooooo---------------------#

  //intersection de la droite avec le premier et le 4eme plan : toutes les longueurs sont en mm
  double theta_deg = 10; //angle en degre
  double theta_rad = theta_deg * TMath::Pi()/180; //angle en radian
		
  double Cible_z = 125.0;
  double dist_L = 100;
  double d = ((Cible_z+dist_L)/cos(theta_rad));
  double d_detect = 14; //distance de l'entree de la boite contenant les detecteurs au premier d√©tecteur
  double moy = 1.25/2; //distance entre 2 detecteurs = 1,25mm, on se place au milieu pour avoir l'equation du plan
  double d_plan4 = 56.75; //7eme detecteur situe a 56.75mm du 1er detecteur  
  
 //intersection avec le 1er plan :
  double Dist1 = d + d_detect + moy;
  double t1 = (Dist1+parFit[0]*sin(theta_rad)-parFit[2]*cos(theta_rad))/(parFit[3]*cos(theta_rad)-parFit[1]*sin(theta_rad));
  pointB1[0] = parFit[0] + parFit[1]*t1;
  cout<<"----------pointB1[0] : "<< pointB1[0]<<endl;
  pointB1[1] = parFit[2] + parFit[3]*t1;
  pointB1[2] = t1;
  
  //intersection avec le 4eme plan :
  double Dist2 = d + d_detect + moy + d_plan4;
  double t2 = (Dist2+parFit[0]*sin(theta_rad)-parFit[2]*cos(theta_rad))/(parFit[3]*cos(theta_rad)-parFit[1]*sin(theta_rad));
  pointB2[0] = parFit[0] + parFit[1]*t2;
  cout<<"----------pointB2[0] : "<< pointB2[0]<<endl;
  pointB2[1] = parFit[2] + parFit[3]*t2;
  pointB2[2] = t2;
 
	}
}*/

void TrackAnalytique::Tracking()
{
  //cout<<"----- In Tracking()"<<endl;
/*  
    if(size==2)
	{
		pointB1[0] = xa[0];
		pointB1[1] = ya[0];
		pointB1[2] = za[0];
		pointB2[0] = xa[1];
		pointB2[1] = ya[1];
		pointB2[2] = za[1];
	}
*/	
	//if(size>=2)
//	{
  //-------------------oooooOOOOO00000OOOOOooooo---------------------#
  //                                            
  //	THE PARAMETRIC LINE EQUATION	  
  //                                       
  //-------------------oooooOOOOO00000OOOOOooooo---------------------#

  // a parameteric line is define from 6 parameters but 4 are independent
  // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
  // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1; 
  // x = p[0] + p[1]*t; 
  // y = p[2] + p[3]*t;
  // z = t;  


  //-------------------oooooOOOOO00000OOOOOooooo---------------------#
  //                                            
  //	EQUATIONS TO SOLVE WITH "LEAST SQUARES" METHOD
  //                                       
  //-------------------oooooOOOOO00000OOOOOooooo---------------------#

  // x = p[0] + p[1]*z; -> equation 1
  // y = p[2] + p[3]*z; -> equation 2
  
 
double xsum, ysum, zsum_x,zsum_y,zzsum_x, zzsum_y, xzsum, yzsum, xsigma2sum, ysigma2sum;

//cout<<"----- Variables created"<<endl;

xsum = 0.0;
ysum = 0.0;
zsum_x = 0.0;
zsum_y = 0.0;
zzsum_x = 0.0;
zzsum_y = 0.0;
xzsum = 0.0;
yzsum = 0.0;
xsigma2sum = 0.0;
ysigma2sum = 0.0;

//cout<<"-----Xa.size : "<<Xa.size()<<endl;
//cout<<"-----Ya.size : "<<Ya.size()<<endl;
//cout<<"-----Za.size : "<<Za.size()<<endl;
//cout<<"-----errorXa.size : "<<errorXa.size()<<endl;
//cout<<"-----errorXa.size : "<<errorYa.size()<<endl;
///cout<<"-----errorXa.size : "<<errorZa.size()<<endl;


for (int i=0;i<Xa.size();i++)
{
  //cout<<"----- Entrance in the loop"<<endl;
  //for equation 1
  zsum_x = zsum_x + Za.at(i)/(errorXa.at(i)*errorXa.at(i));
  //cout<<"----- zsum_x"<<endl;
  xsigma2sum = xsigma2sum + 1/(errorXa.at(i)*errorXa.at(i));
  xsum = xsum + Xa.at(i)/(errorXa.at(i)*errorXa.at(i));   
  zzsum_x = zzsum_x + (Za.at(i)*Za.at(i))/(errorXa.at(i)*errorXa.at(i));
  xzsum = xzsum + Xa.at(i)*Za.at(i)/(errorXa.at(i)*errorXa.at(i));
  
  
  //for equation 2
  zsum_y = zsum_y + Za.at(i)/(errorYa.at(i)*errorYa.at(i)); 
  ysigma2sum = ysigma2sum + 1/(errorYa.at(i)*errorYa.at(i));
  ysum = ysum + Ya.at(i)/(errorYa.at(i)*errorYa.at(i));  
  zzsum_y = zzsum_y + (Za.at(i)*Za.at(i))/(errorYa.at(i)*errorYa.at(i));  
  yzsum = yzsum + Ya.at(i)*Za.at(i)/(errorYa.at(i)*errorYa.at(i));  
  
}

double A1 = zsum_x;
double B1 = xsigma2sum;
double C1 = xsum;
double D1 = zzsum_x;
double E1 = xzsum;

double A2 = zsum_y;
double B2 = ysigma2sum;
double C2 = ysum;
double D2 = zzsum_y;
double E2 = yzsum;

//cout<<setw(15)<<"A1"<<setw(15)<<"B1"<<setw(15)<<"C1"<<setw(15)<<"D1"<<setw(15)<<"E1"<<endl;
//cout<<setw(15)<<A1<<setw(15)<<B1<<setw(15)<<C1<<setw(15)<<D1<<setw(15)<<E1<<endl<<endl;
//cout<<setw(15)<<"A2"<<setw(15)<<"B2"<<setw(15)<<"C2"<<setw(15)<<"D2"<<setw(15)<<"E2"<<endl;
//cout<<setw(15)<<A2<<setw(15)<<B2<<setw(15)<<C2<<setw(15)<<D2<<setw(15)<<E2<<endl<<endl;


parFit[0] = (D1*C1-E1*A1)/(B1*D1-A1*A1);
parFit[1] = (E1*B1-C1*A1)/(B1*D1-A1*A1);
parFit[2] = (D2*C2-E2*A2)/(B2*D2-A2*A2);
parFit[3] = (E2*B2-C2*A2)/(B2*D2-A2*A2);

  //for (int i = 0; i <4; ++i) 
  //{
	//cout<<"----------parFit["<<i<<"] : "<<parFit[i]<<endl;
  //}


  //-------------------oooooOOOOO00000OOOOOooooo---------------------#
  //                                            
  //	FITTED LINE'S INTERSECTION WITH THE FIRST AND THE LAST PLANES
  //                                       
  //-------------------oooooOOOOO00000OOOOOooooo---------------------#

  //intersection de la droite avec le premier et le 4eme plan : toutes les longueurs sont en mm
  double theta_deg = 10; //angle en degre
  double theta_rad = theta_deg * TMath::Pi()/180; //angle en radian
		
  double Cible_z = 125.0;
  double dist_L = 100;
  double d = ((Cible_z+dist_L)/cos(theta_rad));
  double d_detect = 14; //distance de l'entree de la boite contenant les detecteurs au premier d√©tecteur
  double moy = 1.25/2; //distance entre 2 detecteurs = 1,25mm, on se place au milieu pour avoir l'equation du plan
  double d_plan4 = 56.75; //7eme detecteur situe a 56.75mm du 1er detecteur  
  double DetectB_z = 0.025; //demi epaisseur d'un detecteur
  
 //intersection avec le 1er plan :
  double Dist1 = d + d_detect + moy + DetectB_z;
  double t1 = (Dist1+parFit[0]*sin(theta_rad))/(cos(theta_rad)-parFit[1]*sin(theta_rad));
  pointB1[0] = parFit[0] + parFit[1]*t1;
  //cout<<"----------pointB1[0] : "<< pointB1[0]<<endl;
  pointB1[1] = parFit[2] + parFit[3]*t1;
  //cout<<"----------pointB1[1] : "<< pointB1[1]<<endl;
  pointB1[2] = t1;
  //cout<<"----------pointB1[2] : "<< pointB1[2]<<endl;
  
  //intersection avec le 4eme plan :
  double Dist2 = d + d_detect + moy + DetectB_z + d_plan4;
  double t2 = (Dist2+parFit[0]*sin(theta_rad))/(cos(theta_rad)-parFit[1]*sin(theta_rad));
  pointB2[0] = parFit[0] + parFit[1]*t2;
  //cout<<"----------pointB2[0] : "<< pointB2[0]<<endl;
  pointB2[1] = parFit[2] + parFit[3]*t2;
  pointB2[2] = t2;
  
	//}
}


void TrackAnalytique::Line(double t, double *p, double &x, double &y, double &z) 
{ 
  	//cout<<"----------Passage dans line() OK ----------"<<endl;
   // a parameteric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1; 
   x = p[0] + p[1]*t; 
   y = p[2] + p[3]*t;
   z = t; 
} 

TPolyLine3D* TrackAnalytique::GetPolyLine()
{
// draw the fitted line
   double x,y;
   double z0 = -200;
   double z = z0, zmax = 300;
   int n = 100;
   double dt = (zmax - z0) / n;
   int i = 0;
   TPolyLine3D *l = new TPolyLine3D(n);
   while(z < zmax) {
      double t = z0+ dt*i;
      Line(t,parFit,x,y,z);
      l->SetPoint(i,x,y,z);
	  i++;
   }
   return l;
}

double* TrackAnalytique::GetPointB1()
{
  //cout<<"----------Passage dans GetPointB1() OK ----------"<<endl;
  return pointB1;
}

double* TrackAnalytique::GetPointB2()
{ 
 //cout<<"----------Passage dans GetPointB2() OK ----------"<<endl;
 return pointB2;
}

int TrackAnalytique::GetSize()
{
  return size;
}