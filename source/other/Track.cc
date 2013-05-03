#include "Track.hh"



Track::Track()
{
 //first = true;  
 	pointB1 = new double[3];
	pointB2 = new double[3]; 
	pointB1[0] = 0;
	pointB1[1] = 0;
	pointB1[2] = 0;
	pointB2[0] = 1;
	pointB2[1] = 1;
	pointB2[2] = 1;	
	//first = true;
}

Track::Track(vector<double> X, vector<double> Y, vector<double> Z)
//Track::Track(vector<Double_t> X, vector<Double_t> Y, vector<Double_t> Z)
{
  	pointB1 = new double[3];
	pointB2 = new double[3]; 
	pointB1[0] = 0;
	pointB1[1] = 0;
	pointB1[2] = 0;
	pointB2[0] = 1;
	pointB2[1] = 1;
	pointB2[2] = 1;	
  //const Int_t np = (Int_t) X.size();//array's size
  //Double_t* xa = new Double_t [np];
  //Double_t* ya = new Double_t [np];
  //Double_t* za = new Double_t [np];

  size = (Int_t) X.size();
  xa = new Double_t [size];
  ya = new Double_t [size];
  za = new Double_t [size]; 
  
  for (Int_t j=0 ; j<size ; j++)
  {
	xa[j]=X.at(j);
	ya[j]=Y.at(j);
	za[j]=Z.at(j);
  }
  
  cout<<"----------Construction de l'objet track OK----------"<<endl;
  
  //gr = 0;
  //cout<<"----------Initialisation du pointeur gr OK----------"<<endl;

  //first = true; 
  //cout<<"----------Initialisation de first OK----------"<<endl;
}

Track::~Track()
{

	
	cout<<"----------Passage destructeur track OK----------"<<endl;
	//delete [] parFit;
	//delete gr;
	//cout<<"----------delete gr dans destructeur OK----------"<<endl;
}



  //-------------------oooooOOOOO00000OOOOOooooo---------------------#
  //                                            
  //	TOOLS FOR THE LINE 3D FIT		  
  //                                       
  //-------------------oooooOOOOO00000OOOOOooooo---------------------#
  
// define the parameteric line equation 
void Track::line(double t, double *p, double &x, double &y, double &z) 
{ 
  	cout<<"----------Passage dans line() OK ----------"<<endl;
   // a parameteric line is define from 6 parameters but 4 are independent
   // x0,y0,z0,z1,y1,z1 which are the coordinates of two points on the line
   // can choose z0 = 0 if line not parallel to x-y plane and z1 = 1; 
   x = p[0] + p[1]*t; 
   y = p[2] + p[3]*t;
   z = t; 
} 

/*double Track::distance2(double x,double y,double z, double *p) 
{ 
  	cout<<"----------Passage dans distance2() OK ----------"<<endl;
   // distance line point is D= | (xp-x0) cross  ux | 
   // where ux is direction of line and x0 is a point in the line (like t = 0) 
   XYZVector xp(x,y,z); //XYZVector : Class describing a generic displacement vector in 3 dimensions.
						//This class is templated on the type of Coordinate system.
   XYZVector x0(p[0], p[2], 0. ); 
   XYZVector x1(p[0] + p[1], p[2] + p[3], 1. ); 
   XYZVector u = (x1-x0).Unit(); //return unit vector parallel to this->vecteur directeur et unitaire de la droite
   double d2 = ((xp-x0).Cross(u)) .Mag2(); //Cross() : Return vector (cross) product of two displacement vectors,
											//as a vector in the coordinate system of this class.->vector product (!= scalar product)
											//Mag2() : return fX*fX + fY*fY + fZ*fZ, i.e. get magnitude squared
   return d2; 
}*/

// function to be minimized 


void Track::Tracking()
{
  if(size==2)
	{
		pointB1[0] = xa[0];
		pointB1[1] = ya[0];
		pointB1[2] = za[0];
		pointB2[0] = xa[1];
		pointB2[1] = ya[1];
		pointB2[2] = za[1];
		//pointB2[3]={x[1],y[1],z[1]};
	}
	
	if(size>2)
	{
  cout<<"----------Passage dans tracking() OK ----------"<<endl;
    //double e = 0.1;
   Int_t nd = 10000;


//    double xmin = 0; double ymin = 0;
//    double xmax = 10; double ymax = 10;

	TGraph2D * gr = new TGraph2D("position","position des points d'entree des particules secondaires dans les detecteurs",size,xa,ya,za);

   
   // fit the graph now 
   
   TVirtualFitter *min = TVirtualFitter::Fitter(0,4);//static TVirtualFitter*Fitter(TObject* obj, Int_t maxpar = 25)
   min->SetObjectFit(gr);
   min->SetFCN(SumDistance2);//To set the address of the minimization objective function
							 //this function is called by CINT instead of the function above
  Double_t arglist[10];
  arglist[0] = 3;
  //min->ExecuteCommand("SET PRINT",arglist,1);
  
    
   double pStart[4] = {1,1,1,1};
   min->SetParameter(0,"x0",pStart[0],0.000001,0,0); // Last 3 parameters : step size, min, max (min=max=0 => no limits)
   min->SetParameter(1,"Ax",pStart[1],0.000001,0,0);
   min->SetParameter(2,"y0",pStart[2],0.000001,0,0);
   min->SetParameter(3,"Ay",pStart[3],0.000001,0,0);
   
  arglist[0] = 1000; // number of function calls 
  arglist[1] = 0.0000001; // tolerance 
  
  min->ExecuteCommand("MIGRAD",arglist,2); // Minimization method
  
    //if (minos) min->ExecuteCommand("MINOS",arglist,0);
   int nvpar,nparx; 
   double amin,edm, errdef;
   min->GetStats(amin,edm,errdef,nvpar,nparx);
   min->PrintResults(1,amin);
   
     // get fit parameters
  double parFit[4];
  for (int i = 0; i <4; ++i) 
  {
	parFit[i] = min->GetParameter(i);
	cout<<"----------parFit["<<i<<"] = min->GetParameter("<<i<<") OK----------"<<endl;
  }
  

  
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
  double D1 = d + d_detect + moy;
  double t1 = (D1+parFit[0]*sin(theta_rad)-parFit[2]*cos(theta_rad))/(parFit[3]*cos(theta_rad)-parFit[1]*sin(theta_rad));
  pointB1[0] = parFit[0] + parFit[1]*t1;
  cout<<"----------pointB1[0] OK----------"<<endl;
  pointB1[1] = parFit[2] + parFit[3]*t1;
  pointB1[2] = t1;
  
  //intersection avec le 4eme plan :
  double D2 = d + d_detect + moy + d_plan4;
  double t2 = (D2+parFit[0]*sin(theta_rad)-parFit[2]*cos(theta_rad))/(parFit[3]*cos(theta_rad)-parFit[1]*sin(theta_rad));
  pointB2[0] = parFit[0] + parFit[1]*t2;
  cout<<"----------pointB2[0] OK----------"<<endl;
  pointB2[1] = parFit[2] + parFit[3]*t2;
  pointB2[2] = t2;
   
  delete gr;
  cout<<"----------delete gr dans tracking OK----------"<<endl;
	}
}




 

double* Track::GetPointB1()
{
  cout<<"----------Passage dans GetPointB1() OK ----------"<<endl;
  return pointB1;
}

double* Track::GetPointB2()
{ 
 cout<<"----------Passage dans GetPointB2() OK ----------"<<endl;
 return pointB2;
}


void SumDistance2(int &, double *, double & sum, double * par, int )
{   
//  Track MyTrack;
  bool first = true; 
  // the TGraph must be a global variable
  TGraph2D * gr = dynamic_cast<TGraph2D*>( (TVirtualFitter::GetFitter())->GetObjectFit() );//GetFitter : return the current Fitter
																							//fObjectFit : pointer to object being fitted
  assert(gr != 0);
  double * x = gr->GetX();
  double * y = gr->GetY();
  double * z = gr->GetZ();
  int npoints = gr->GetN();
  sum = 0;
  for (int i  = 0; i < npoints; ++i) 
  { 
//	  double d = MyTrack.distance2(x[i],y[i],z[i],par); 
	  double d = distance2(x[i],y[i],z[i],par); 
	  sum += d;
	  #ifdef DEBUG
	  if (first) std::cout << "point " << i << "\t" 
						  << x[i] << "\t" 
						  << y[i] << "\t" 
						  << z[i] << "\t" 
						  << std::sqrt(d) << std::endl; 
	  #endif
  }
  //if (first) std::cout << "Total sum2 = " << sum << std::endl;
  first = false;
}

double distance2(double x,double y,double z, double *p) 
{ 
  	//cout<<"----------Passage dans distance2() OK ----------"<<endl;
   // distance line point is D= | (xp-x0) cross  ux | 
   // where ux is direction of line and x0 is a point in the line (like t = 0) 
   XYZVector xp(x,y,z); //XYZVector : Class describing a generic displacement vector in 3 dimensions.
						//This class is templated on the type of Coordinate system.
   XYZVector x0(p[0], p[2], 0. ); 
   XYZVector x1(p[0] + p[1], p[2] + p[3], 1. ); 
   XYZVector u = (x1-x0).Unit(); //return unit vector parallel to this->vecteur directeur et unitaire de la droite
   double d2 = ((xp-x0).Cross(u)) .Mag2(); //Cross() : Return vector (cross) product of two displacement vectors,
											//as a vector in the coordinate system of this class.->vector product (!= scalar product)
											//Mag2() : return fX*fX + fY*fY + fZ*fZ, i.e. get magnitude squared
   return d2; 
}

