  /// \file Calorimeter.hh
/*
 *
 * Calorimeter.hh header template generated by fclass
 * Creation date : mer. juin 12 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef CALORIMETER_HH
#define CALORIMETER_HH

#include "Detector/Detector.hh"

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector>
#include <algorithm>


#include "Objects/CaloHit.hh"

namespace baboon {


	enum CalorimeterType {

		fHCAL,
		fECAL,
		fLumiCAL,
		fBCAL,
		fLHCAL,
		fCustomCAL
	};
	/*!
	 * Class Calorimeter
	 * Inherits from base class Detector
	 */

	class Calorimeter : public Detector {

		public:

			/*!
			 *
			 * Default Constructor
			 *
			 */
			Calorimeter( const CalorimeterType &calType );

			/*!
			 *
			 * Default Destructor
			 *
			 */
			virtual ~Calorimeter();

			/*!
			 *
			 *
			 *
			 */
			virtual Return AddCaloHit( CaloHit *caloHit );

			/*!
			 *
			 *
			 *
			 */
			virtual Return RemoveCaloHit( CaloHit *caloHit );



		protected:

			CalorimeterType calorimeterType;
			CaloHitCollection *caloHitCollection;

		public:

			inline CalorimeterType GetType()
				{ return calorimeterType; }




	};  // class

}  // namespace 

#endif  //  CALORIMETER_HH
