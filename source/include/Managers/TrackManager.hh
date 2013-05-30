  /// \file TrackManager.hh
/*
 *
 * TrackManager.hh header template generated by fclass
 * Creation date : jeu. mai 9 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef TRACKMANAGER_HH
#define TRACKMANAGER_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector>
#include <algorithm>

#include "Reconstruction/Objects/Track.hh"
#include "Utilities/ReturnValues.hh"
#include "Utilities/Globals.hh"


namespace baboon {

	/*
	 * Class TrackManager
	 */

	class TrackManager {

		public:

			static TrackManager *GetInstance();

			static void Kill();

			Return AddTrack( Track *track );

			Return RemoveTrack( Track *track );

			bool TrackContainsHit( Track *track , Hit *hit );

			Return ClearAllContent();

			Return FindTrackContainingHit( Hit *hit , Track *trackToFind );


		protected:

			/*! Default Constructor */
			TrackManager();
			/*! Default Destructor */
			virtual ~TrackManager();


			static TrackManager *instance;
			TrackCollection *trackCollection;

		public:

			inline TrackCollection *GetTrackCollection()
				{ return trackCollection; }



	};  // class

}  // namespace 

#endif  //  TRACKMANAGER_HH
