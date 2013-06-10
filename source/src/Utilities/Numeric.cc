/*
 *
 * Numeric.cc cpp file template generated by fclass
 * Creation date : Fri Mar 15 14:10:39 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#include "Utilities/Numeric.hh"



using namespace std ;



namespace baboon {


	double square(double x,double y) {
		return x += y*y;
	}


	std::pair<double,double> linearRegression( vector<double>& x , vector<double>& y ) {

		double avX = accumulate( x.begin() , x.end() , 0.0 ) *1.0 / x.size();
		double avY = accumulate( y.begin() , y.end() , 0.0 ) *1.0 / y.size();
		double avSquareX = accumulate( x.begin() , x.end() , 0.0 , square);

		double squaredSum = 0;

		for(unsigned int i=0 ; i<x.size() ; i++) {
			squaredSum += x.at(i)*y.at(i);
		}

		double a = squaredSum - avX*avY;
		a /= (avSquareX - avX*avX);

		double b = avY - a*avX;
		return make_pair(a,b);
	}

	double gauss2D(double *x, double *par) {

	   double z1 = double((x[0]-par[1])/par[2]);
	   double z2 = double((x[1]-par[3])/par[4]);

	   return par[0]*exp(-0.5*(z1*z1+z2*z2));

	}

	double doubleGauss2D(double *x, double *par) {

	   return gauss2D(x,&par[0]) + gauss2D(x,&par[5]);
	}



	double Mean( const std::vector<double> &vec ) {

		return accumulate( vec.begin() , vec.end() , 0.0 ) *1.0 / vec.size();
	}



	double RMS( const std::vector<double> &vec ) {

		double mean = Mean( vec );
		double meanSquare = accumulate( vec.begin() , vec.end() , 0.0 , square);
		return sqrt( (meanSquare - mean*mean ) / vec.size() );
	}



	double RMS90( const std::vector<double> &vec ) {

		vector<double> distanceFromMean;
		int erasePosition = int( double(vec.size()) / 10.0 );
		double mean = Mean( vec );

		for( unsigned int i=0 ; i<vec.size() ; i++ )
			distanceFromMean.push_back( abs(mean - vec.at(i))  );

		std::sort( distanceFromMean.begin() , distanceFromMean.end() );
		distanceFromMean.erase( distanceFromMean.begin() , distanceFromMean.begin() + erasePosition );

		return RMS( distanceFromMean );
	}

	double RecoveryProbabilityWithinSigma( const std::vector<double> &vec , const int &nbOfSigmas ) {

		double mean = Mean( vec );
		double sigma = RMS( vec );
		int entriesWithinSigmas = 0;
		for( unsigned int i=0 ; i<vec.size() ; i++ )
			if( abs(vec.at(i) - mean ) < sigma*nbOfSigmas )
				entriesWithinSigmas++;

		return double( vec.size() ) / double( entriesWithinSigmas ) ;
	}


}
