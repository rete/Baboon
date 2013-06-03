  /// \file ConeBeginningAlgorithm.hh
/*
 *
 * ConeBeginningAlgorithm.hh header template generated by fclass
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


#ifndef CONEBEGINNINGALGORITHM_HH
#define CONEBEGINNINGALGORITHM_HH

#include "Algorithm/AbstractAlgorithm.hh"
#include "Exception.hh"

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 
#include <algorithm>
#include <numeric>

#include "Geometry/Cone.hh"
#include "Objects/HitCollection.hh"
#include "Managers/HitManager.hh"
#include "Managers/CoreManager.hh"
#include "Managers/ShowerManager.hh"
#include "Managers/ClusteringManager.hh"
#include "Objects/Shower.hh"
#include "Utilities/ReturnValues.hh"
#include "Utilities/Globals.hh"

#include "Geometry/SDHCALPrototype.hh"
#include "TGeoCone.h"

namespace baboon {

	/*!
	 * @brief Class ConeBeginningAlgorithm
	 * Inherits from base class AbstractAlgorithm
	 */

	class ConeBeginningAlgorithm : public AbstractAlgorithm {

		public:

			/*!
			 *
			 * @brief Default Constructor
			 *
			 */
			ConeBeginningAlgorithm();

			/*!
			 *
			 * @brief Default Destructor
			 *
			 */
			virtual ~ConeBeginningAlgorithm();


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
			 * @brief Find the existing pads inside a given cone. Find also the touched pads.
			 *
			 */
			Return FindPadsInCone( Cone *cone , std::vector<IntVector> &existingPadsInCone, std::vector<IntVector> &touchedPadsInCone );


			int NbOfCoreHitsInCone( Cone *cone , Core *core );

			HitCollection *hitCollection;
			std::vector<TGeoNode*> nodeList;

			// algorithm parameters
			int xLoopHalfRange;
			int yLoopHalfRange;
			int zLoopHalfRange;
			double coneMinimumOpeningAngle;
			double coneMaximumOpeningAngle;
			double angleStep;
			int coneLength;
			double coneDirectedVectorHalfRange;
			double coneDirectedVectorStep;
			double conePadsLowerLimitCut;
			int coreSizeLowerLimit;
			double coneBackwardDistance;
			double coneStepDistance;
			int clusterSizeAtPeakPositionUpperLimit;


	};  // class


}  // namespace 


#endif  //  CONEBEGINNINGALGORITHM_HH
