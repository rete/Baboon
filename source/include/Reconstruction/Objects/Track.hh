  /// \file Track.hh
/*
 *
 * Track.hh header template generated by fclass
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


#ifndef TRACK_HH
#define TRACK_HH

#include "Reconstruction/Objects/RecoObject.hh"

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector>
#include <algorithm>


#include "Utilities/ReturnValues.hh"
#include "Objects/Hit.hh"
#include "Geometry/ThreeVector.hh"

namespace baboon {




	typedef std::pair<Hit*,Hit*> TrackExtremities;


	/*
	 * Class Track
	 * Inherits from base class RecoObject
	 */

	class Track : public RecoObject {

		public:

			/*!
			 * @brief : Default Constructor
			 */
			Track();

			/*!
			 * @brief : Default Destructor
			 */
			virtual ~Track();

			/*!
			 *
			 * @brief : Return the hit collection which define the track segment
			 *
			 */
			inline HitCollection *GetHitCollection()
				{ return hitCollection; }


			TrackExtremities GetExtremities();

			double GetTrackLength();

			void SetTrackAjustementVector( const ThreeVector& trackVec)
				{ trackAjustementVector = trackVec; }

			std::vector<ThreeVector> GetPositions() const;

			Return AddHit( Hit *hit );

			Return RemoveHit( Hit *hit );

			inline unsigned int Size()
				{ return hitCollection->size(); }

		protected:

			HitCollection *hitCollection;
			TrackExtremities extremities;
			ThreeVector trackAjustementVector;


	};  // class

	typedef std::vector<Track*> TrackCollection;

}  // namespace 

#endif  //  TRACK_HH
