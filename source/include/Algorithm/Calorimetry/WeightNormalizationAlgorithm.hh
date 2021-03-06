  /// \file WeightNormalizationAlgorithm.hh
/*
 *
 * WeightNormalizationAlgorithm.hh header template generated by fclass
 * Creation date : ven. juil. 19 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef WEIGHTNORMALIZATIONALGORITHM_HH
#define WEIGHTNORMALIZATIONALGORITHM_HH

#include "Algorithm/AbstractAlgorithm.hh"

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 

#include "Managers/ShowerManager.hh"
#include "Objects/Shower.hh"
#include "Utilities/CaloHitHelper.hh"

namespace baboon {

	/*!
	 *
	 * @brief  Class WeightNormalizationAlgorithm
	 * 		   Inherits from base class AbstractAlgorithm
	 *
	 */

	class WeightNormalizationAlgorithm : public AbstractAlgorithm {

		public:

			struct HitShowerWeight {
				Shower *shower;
				CaloHit *caloHit;
				double weight;
			};

			typedef std::vector<HitShowerWeight*> HitShowerWeightCollection;

			/*!
			 *
			 * Default Constructor
			 *
			 */
			WeightNormalizationAlgorithm();

			/*!
			 *
			 * Default Destructor
			 *
			 */
			virtual ~WeightNormalizationAlgorithm();

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


	};  // class

}  // namespace 

#endif  //  WEIGHTNORMALIZATIONALGORITHM_HH
