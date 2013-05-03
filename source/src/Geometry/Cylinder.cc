/*
 *
 * Cylinder.cc cpp file template generated by fclass
 * Creation date : Wed Mar 27 10:50:03 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#include "Geometry/Cylinder.hh"


using namespace std ;



namespace baboon {


	Cylinder::Cylinder(ThreeVector x1 , ThreeVector x2 , double r) {

		centerBase1 = x1;
		centerBase2 = x2;
		radius = r;
		directedVector = x2 - x1;
		directedVector.setMag(1);
		length = (x2 - x1).mag();
		ThreeVector temp = x2 - x1;
		temp.setMag( temp.mag() / 2.0 );
		center = x1 + temp;
	}

	Cylinder::Cylinder(ThreeVector c , ThreeVector vecDir , double l , double r) {

		center = c;
		directedVector = vecDir;
		length = l;
		radius = r;
		ThreeVector temp = vecDir;
		temp.setMag( length / 2.0 );
		centerBase1 = c - temp;
		centerBase2 = c + temp;

	}

	Cylinder::~Cylinder() {}

	void Cylinder::SetLength( double l ) {

		if( l > 0 ) {
			ThreeVector temp = directedVector * ( l/2.0 - length/2.0 );
			centerBase1 = centerBase1 - temp;
			centerBase2 = centerBase2 + temp;
			length = l;
		}
		else {
			length = 0;
			centerBase1 = centerBase2 = center;
		}
	}

	void Cylinder::SetRadius( double r ) {

		if( r > 0 ) radius = r;
		else r = 0;

	}

	void Cylinder::SetCenterBase1( const ThreeVector& vec1 ) {

		centerBase1 = vec1;
		ThreeVector temp = centerBase2 - centerBase1;
		length = temp.mag();
		directedVector = temp;
		directedVector.setMag(1);
		temp.setMag( length/2.0 );
		center = centerBase1 + temp;
	}

	void Cylinder::SetCenterBase2( const ThreeVector& vec2 ) {

		centerBase2 = vec2;
		ThreeVector temp = centerBase2 - centerBase1;
		length = temp.mag();
		directedVector = temp;
		directedVector.setMag(1);
		temp.setMag( length/2.0 );
		center = centerBase1 + temp;
	}

	void Cylinder::SetCenter( const ThreeVector& c ) {

		ThreeVector temp = c - center;
		centerBase1 = centerBase1 + temp;
		centerBase2 = centerBase2 + temp;
		center = c;

	}


	void Cylinder::SetDirectedVector( const ThreeVector& dirVec) {

		directedVector = dirVec;
		ThreeVector temp = dirVec;
		temp.setMag( length / 2.0 );
		centerBase1 = center - temp;
		centerBase2 = center + temp;
	}


	ThreeVector Cylinder::GetProjectionOnAxe( const ThreeVector& v ) {

		ThreeVector u = centerBase2 - centerBase1;
		u.setMag( (v-centerBase1).mag() * cos( u.angle( v-centerBase1 ) ) );
		return (centerBase1 + u);
	}


	bool Cylinder::Contains( const ThreeVector& v ) {

		ThreeVector projection = this->GetProjectionOnAxe(v);
		if( (projection - v).mag() > radius ) return false;
		else if( (projection - centerBase1).mag() > length || (projection - centerBase2).mag() > length ) return false;
		return true;
	}

}
