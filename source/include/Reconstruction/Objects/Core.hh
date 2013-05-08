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

#include "Reconstruction/Objects/RecoObject.hh"

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector>
#include <algorithm>

#include "Objects/Hit.hh"
#include "Managers/HitManager.hh"
#include "Utilities/ReturnValues.hh"

namespace baboon {


	class Core;


	typedef std::vector<Core*> CoreCollection;

	/*!
	 *
	 * @brief Class Core
	 * Inherits from base class RecoObject
	 *
	 */
	class Core : public RecoObject {

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
			 * @brief Add a hit to the hit collection of the core
			 *
			 */
			Return AddHit( Hit *hit );

			/*!
			 *
			 * @brief Remove a hit from the hit collection of the core
			 *
			 */
			Return RemoveHit( Hit *hit );

		protected:

			HitCollection *hitCollection;


	};  // class

}  // namespace 

#endif  //  CORE_HH
