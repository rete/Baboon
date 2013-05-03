  /// \file Cone.cc
/*
 *
 * Cone.cc source template generated by fclass
 * Creation date : ven. avr. 19 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Geometry/Cone.hh"

namespace sdhcal {

	Cone::Cone( const ThreeVector &peakPos , double thet , const ThreeVector &u , double r ) {

		if( r < 0 ) radius = 0;
		else radius = r;
		if( thet < 0 ) theta = 0;
		else theta = thet;
		peakPosition = peakPos;
		direcVector = u;
	}


	Cone::Cone( const ThreeVector &peakPos , const ThreeVector &v , double r ) {

		peakPosition = peakPos;
		direcVector = v; direcVector.setMag(1);
		if( r < 0 ) radius = 0;
		else radius = r;
		theta = atan(r/v.mag());

	}

	Cone::~Cone() {

	}



	void Cone::SetPeakPosition( const ThreeVector &p ) {
		peakPosition = p;
	}


	void Cone::SetTheta( double t ) {
		if( t < 0 ) theta = 0;
		else theta = t;
	}


	void Cone::SetDirectedVector( const ThreeVector &u ) {
		direcVector = u;
	}


	void Cone::SetRadius( double r ) {
		if( r < 0 ) radius = 0;
		else radius = r;
	}


	bool Cone::Contains( const ThreeVector &v ) {

		ThreeVector posInCone = v-peakPosition;
		if( posInCone.angle(direcVector) > theta ) return false;
		else {
			double h = radius / tan(theta);
			if( direcVector.angle(direcVector*h - posInCone) > M_PI/2.0 ) return false;
		}

		return true;
	}


}  // namespace 

