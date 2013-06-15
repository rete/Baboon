  /// \file Detector.hh
/*
 *
 * Detector.hh header template generated by fclass
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


#ifndef DETECTOR_HH
#define DETECTOR_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector>

#include "Utilities/ReturnValues.hh"

namespace baboon {

	/*!
	 * Class Detector
	 */

	class Detector {

		public:

		   /*!
			*
			* Default Constructor
			*
			*/
			Detector( const std::string &detName = "" );

		   /*!
			*
			* Default Destructor
			*
			*/
			virtual ~Detector();

			virtual Return Build() = 0;

		protected:

			std::string detectorName;



};  // class 

}  // namespace 

#endif  //  DETECTOR_HH