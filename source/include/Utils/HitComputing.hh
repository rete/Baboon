  /// \file HitComputing.hh
/*
 *
 * HitComputing.hh header template generated by fclass
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


#ifndef HITCOMPUTING_HH
#define HITCOMPUTING_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 

#include "EVENT/LCCollection.h"

#include "Objects/Hit.hh"
#include "Geometry/ThreeVector.hh"


namespace baboon {

	/*! Return the closest hit from position v ( in x,y,z coordinates ) */
	Hit *GetClosestHit( const ThreeVector &v , HitCollection *hitCol );

	/*! Return the RMS dispersion of the hit collection */
	double GetRMSDispersion( HitCollection *hitCol );

	/*! Return the center of gravity of the hit collection pondered by thresholds */
	ThreeVector GetCenterOfGravity( HitCollection * hitCol );

	/*! Return the center of gravity of the hit collection pondered by thresholds */
	ThreeVector GetCenterOfGravity( LCCollection *lcCollection );


}  // namespace 

#endif  //  HITCOMPUTING_HH
