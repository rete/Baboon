  /// \file Core.hh
/*
 *
 * Core.hh header template generated by fclass
 * Creation date : lun. avr. 29 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef CORE_HH
#define CORE_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector>
#include <algorithm>

#include "Objects/HitCompositeObject.hh"
#include "Utilities/ReturnValues.hh"
#include "Geometry/Cone.hh"

namespace baboon {


	class Core;


	typedef std::vector<Core*> CoreCollection;

	/*!
	 *
	 * @brief Class Core
	 *
	 */
	class Core : public HitCompositeObject , public TypedObject {

		public:

			/*!
			 *
			 * @brief Default Constructor
			 *
			 */
			Core();

			/*!
			 *
			 * @brief Default Destructor
			 *
			 */
			virtual ~Core();

			/*!
			 *
			 *
			 *
			 */
			Return SetBuildConcentration( const double c );

			/*!
			 *
			 *
			 *
			 */
			inline const double GetBuildConcentration()
				{ return buildConcentration; }

		protected:

			double buildConcentration;




	};  // class

}  // namespace 

#endif  //  CORE_HH
