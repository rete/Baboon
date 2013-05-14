  /// \file CoreManager.hh
/*
 *
 * CoreManager.hh header template generated by fclass
 * Creation date : ven. mai 10 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef COREMANAGER_HH
#define COREMANAGER_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 

#include "Reconstruction/Objects/Core.hh"

namespace baboon {

	/*
	 * Class CoreManager
	 */

	class CoreManager {

		public:

			/*!
			*
			* Return a single instance of this class. Singleton pattern
			*
			*/

			static CoreManager *GetInstance();

			/*!
			*
			* Kill the instance carefully. Singleton pattern
			*
			*/
			static void Kill();

			Return AddCore( Core *core );

			Return RemoveCore( Core *core );

			Return ClearAllContent();

			bool CoreContainsHit( Core *core , Hit *hit );

		protected:

			/*!
			*
			* Unique instance of this class. Singleton pattern
			*
			*/
			static CoreManager* instance;

			/*!
			*
			* Protected constructor. Singleton pattern
			*
			*/
			CoreManager();

			/*!
			*
			* Protected constructor. Singleton pattern
			*
			*/
			virtual ~CoreManager();

			CoreCollection *coreCollection;

	};  // class

}  // namespace 

#endif  //  COREMANAGER_HH