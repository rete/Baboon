  /// \file Cone.hh
/*
 *
 * Cone.hh header template generated by fclass
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


#ifndef CONE_HH
#define CONE_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 

#include "Geometry/ThreeVector.hh"


namespace sdhcal {

	/*
	 * Class Cone
	 */

	class Cone {

		public:

			/*!  First constructor version */
			Cone( const ThreeVector &peakPos , double thet , const ThreeVector &u , double r );

			/*!  Second constructor version */
			Cone( const ThreeVector &peakPos , const ThreeVector &v , double r );

			/*! Default Destructor */
			virtual ~Cone();

			/*!  */
			void SetPeakPosition( const ThreeVector & );

			/*!  */
			void SetTheta( double t );

			/*!  */
			void SetDirectedVector( const ThreeVector & );

			/*!  */
			void SetRadius( double r );

			/*!  */
			inline ThreeVector GetPeakPosition()
				{ return peakPosition; }

			/*!  */
			inline double GetTheta()
				{ return theta; }

			/*!  */
			inline ThreeVector GetDirectedVector()
				{ return direcVector; }

			/*!  */
			inline double GetRadius()
				{ return radius; }

			/*! */
			bool Contains( const ThreeVector &v );


		protected:

			ThreeVector peakPosition;

			ThreeVector direcVector;

			double radius;

			double theta;



	};  // class


}  // namespace 

#endif  //  CONE_HH
