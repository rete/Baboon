  /// \file TrackToShowerAssociationAlgorithm.hh
/*
 *
 * TrackToShowerAssociationAlgorithm.hh header template generated by fclass
 * Creation date : lun. mai 27 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : remi
 */


#ifndef TRACKTOSHOWERASSOCIATIONALGORITHM_HH
#define TRACKTOSHOWERASSOCIATIONALGORITHM_HH

#include "Algorithm/AbstractAlgorithm.hh"

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector>


#include "Managers/ShowerManager.hh"
#include "Managers/TrackManager.hh"
#include "Reconstruction/Objects/Track.hh"
#include "Reconstruction/Objects/Shower.hh"
#include "Reconstruction/Linear3DFit.hh"
#include "Geometry/Cylinder.hh"
#include "Utilities/Internal.hh"

//#include "Geometry/SDHCALPrototype.hh"
//#include "TGeoCyl.h"

namespace baboon {

	/*
	 * Class TrackToShowerAssociationAlgorithm
	 * Inherits from base class AbstractAlgorithm
	 */

	class TrackToShowerAssociationAlgorithm : public AbstractAlgorithm {

		public:

			/*!
			 *
			 * Default Constructor
			 *
			 */
			TrackToShowerAssociationAlgorithm();

			/*!
			 *
			 * Default Destructor
			 *
			 */
			virtual ~TrackToShowerAssociationAlgorithm();

		protected:

			/*!
			 *
			 * @brief Initialize the algorithm, i.e by initializing specific variables
			 *
			 */
			virtual void Init();


			/*!
			 *
			 * @brief Execute the algorithm
			 *
			 */
			virtual void Execute();


			/*!
			 *
			 * @brief Finalize the algorithm
			 *
			 */
			virtual void End();


			/*!
			 *
			 * @brief Allow to check if everything is well set in the algorithm before starting it
			 *
			 */
			virtual Return CheckConsistency();

			Return PredominantShowerInCylinder( Cylinder *cyl , Shower *predominantShower );

			Return MergeTrackInShower( Shower *ShowerToEnlarge , Track *trackToAssociate );

			double Chi2LimitForTrack;
			double cylinderLength;
			double cylinderRadius;


	};  // class

}  // namespace 

#endif  //  TRACKTOSHOWERASSOCIATIONALGORITHM_HH
