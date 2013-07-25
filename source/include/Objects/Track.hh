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

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector>
#include <algorithm>


#include "Utilities/ReturnValues.hh"
#include "Objects/HitCompositeObject.hh"
#include "Objects/HitCollection.hh"
#include "Geometry/ThreeVector.hh"

namespace baboon {




	typedef std::pair<Hit*,Hit*> TrackExtremities;


	/*
	 * Class Track
	 */

	class Track : public HitCompositeObject , public TypedObject {

		public:

			/*!
			 * @brief : Default Constructor
			 */
			Track();

			/*!
			 * @brief : Default Destructor
			 */
			virtual ~Track();

			TrackExtremities GetExtremities();

			void SetTrackAjustementVector( const ThreeVector& trackVec)
				{ trackAjustementVector = trackVec; }

			std::vector<ThreeVector> GetPositions() const;

			Return SortHits();

			std::vector<ThreeVector> GetIJKs() const;

		protected:

			TrackExtremities extremities;
			ThreeVector trackAjustementVector;


	};  // class

	typedef std::vector<Track*> TrackCollection;

}  // namespace 

#endif  //  TRACK_HH
