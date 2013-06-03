  /// \file CoreToCoreAssociationAlgorithm.hh
/*
 *
 * CoreToCoreAssociationAlgorithm.hh header template generated by fclass
 * Creation date : ven. mai 24 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef CORETOCOREASSOCIATIONALGORITHM_HH
#define CORETOCOREASSOCIATIONALGORITHM_HH

#include "Algorithm/AbstractAlgorithm.hh"

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 

#include "Managers/CoreManager.hh"
#include "Managers/ShowerManager.hh"
#include "Objects/Core.hh"
#include "Objects/Track.hh"
#include "Objects/Shower.hh"


namespace baboon {

	/*
	 * Class CoreToCoreAssociationAlgorithm
	 * Inherits from base class AbstractAlgorithm
	 */

	class CoreToCoreAssociationAlgorithm : public AbstractAlgorithm {

		public:
			/*!
			 *
			 * Default Constructor
			 *
			 */
			CoreToCoreAssociationAlgorithm();

			/*!
			 *
			 * Default Destructor
			 *
			 */
			virtual ~CoreToCoreAssociationAlgorithm();

		protected:

			/*!
			 *
			 * @brief Initialize the algorithm, i.e by initializing specific variables
			 *
			 */
			virtual Return Init();


			/*!
			 *
			 * @brief Execute the algorithm
			 *
			 */
			virtual Return Execute();


			/*!
			 *
			 * @brief Finalize the algorithm
			 *
			 */
			virtual Return End();


			/*!
			 *
			 * @brief Allow to check if everything is well set in the algorithm before starting it
			 *
			 */
			virtual Return CheckConsistency();

			/*!
			 *
			 * @brief Return true if the cores can be associated in the algorithm sense.
			 *
			 */
			bool CanAssociateCores( Core *core1 , Core *core2 );


	};  // class

}  // namespace 

#endif  //  CORETOCOREASSOCIATIONALGORITHM_HH