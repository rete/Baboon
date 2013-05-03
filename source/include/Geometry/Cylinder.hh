/*
 *
 * Cylinder.hh header template generated by fclass
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


#ifndef CYLINDER_HH
#define CYLINDER_HH

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>


#include "Geometry/ThreeVector.hh"


namespace sdhcal {

	/*!
	 * Class Cylinder.
	 */

	class Cylinder {

		protected :

			double radius;

			double length;

			ThreeVector centerBase1;

			ThreeVector centerBase2;

			ThreeVector directedVector;

			ThreeVector center;


		public :
			// Defined by the position of the base center 1 and 2 and the radius.
			Cylinder(ThreeVector x1 , ThreeVector x2 , double r);

			// Defined by the position in the center, the directed vector of the main axe, the length and the radius.
			Cylinder(ThreeVector c , ThreeVector vecDir , double l , double r);

			/*! Default Destructor */
			virtual ~Cylinder();

			void SetLength( double l );

			void SetRadius( double r );

			void SetCenterBase1( const ThreeVector& );

			void SetCenterBase2( const ThreeVector& );

			void SetCenter( const ThreeVector& );

			void SetDirectedVector( const ThreeVector& );

			inline double GetLength()
				{ return length; }

			inline double GetRadius()
				{ return radius; }

			inline ThreeVector GetCenterBase1()
				{ return centerBase1; }

			inline ThreeVector GetCenterBase2()
				{ return centerBase2; }

			inline ThreeVector GetCenter()
				{ return center; }

			inline ThreeVector GetDirectedVector()
				{ return directedVector; }

			ThreeVector GetProjectionOnAxe( const ThreeVector& );

			bool Contains( const ThreeVector& v );


	};

}

#endif  // CYLINDER_HH
