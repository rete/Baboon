  /// \file PrincipalComponentAnalysis.cc
/*
 *
 * PrincipalComponentAnalysis.cc source template generated by fclass
 * Creation date : dim. avr. 21 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Algorithm/PrincipalComponentAnalysis.hh"

using namespace std;

namespace baboon {

	PrincipalComponentAnalysis::PrincipalComponentAnalysis()
		: AbstractAlgorithm("PrincipalComponentAnalysis") {

		needData = false;

	}

	PrincipalComponentAnalysis::~PrincipalComponentAnalysis() {

	}


	Return PrincipalComponentAnalysis::Init() {

		eigenValues.Clear();
		eigenVectors.Clear();
		return BABOON_SUCCESS();
	}


	Return PrincipalComponentAnalysis::CheckConsistency() {

		if( hitCollection == 0 )
			return BABOON_ERROR("Principal Component Analysis bad init. Please check your inputs!");
		return BABOON_SUCCESS();
	}


	Return PrincipalComponentAnalysis::Execute() {

		double *x = new double[ hitCollection->size() ];
		double *y = new double[ hitCollection->size() ];

		for( unsigned int i=0 ; i<hitCollection->size() ; i++ ) {

			x[i] = hitCollection->at(i)->GetPosition().x();
			y[i] = hitCollection->at(i)->GetPosition().y();

		}

		double xMean = 0;
		double yMean = 0;

		TMatrixD covMatrix(2,2);

		covMatrix(0,0) = 0;
		covMatrix(0,1) = 0;
		covMatrix(1,0) = 0;
		covMatrix(1,1) = 0;

		for( int i=0 ; i <hitCollection->size() ; i++ ) {
			xMean += x[i] / hitCollection->size();
			yMean += y[i] / hitCollection->size();
		}

		for( int i=0 ; i <hitCollection->size() ; i++ ) {

			x[i] -= xMean;
			y[i] -= yMean;
			covMatrix(0,0) += ( x[i]*x[i] ) / (hitCollection->size()-1);
			covMatrix(0,1) += ( x[i]*y[i] ) / (hitCollection->size()-1);
			covMatrix(1,0) += ( y[i]*x[i] ) / (hitCollection->size()-1);
			covMatrix(1,1) += ( y[i]*y[i] ) / (hitCollection->size()-1);

		}

		TMatrixDEigen matrixEigen(covMatrix);

		eigenValues.ResizeTo(2);
		eigenVectors.ResizeTo(2,2);
		eigenValues = matrixEigen.GetEigenValuesRe();
		eigenVectors = matrixEigen.GetEigenVectors();

		delete x;
		delete y;

		return BABOON_SUCCESS();
	}


	Return PrincipalComponentAnalysis::End() {

		return BABOON_SUCCESS();
	}

}  // namespace 













