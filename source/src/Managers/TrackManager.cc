  /// \file TrackManager.cc
/*
 *
 * TrackManager.cc source template generated by fclass
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


#include "Managers/TrackManager.hh"

using namespace std;

namespace baboon {

	TrackManager *TrackManager::instance = 0;


	TrackManager::TrackManager() {

		trackCollection = new TrackCollection();
	}

	TrackManager::~TrackManager() {}

	TrackManager *TrackManager::GetInstance() {

		if( instance ==0 )
			instance = new TrackManager();
		return instance;
	}


	void TrackManager::Kill() {

		if( instance != 0 ) {
			delete instance;
			instance = 0;
		}
	}


	Return TrackManager::AddTrack( Track *track ) {

		if( trackCollection->empty() ) {
			trackCollection->push_back( track );
			return S_OK();
		}
		TrackCollection::iterator trackIt = std::find( trackCollection->begin() , trackCollection->end() , track );

		if( trackIt == trackCollection->end() ) {
			trackCollection->push_back( track );
			return S_OK();
		}

		return S_ERROR("While adding track. Track already registered by the Track Manager!");
	}



	Return TrackManager::RemoveTrack( Track *track ) {

		if( trackCollection->empty() )
			return S_ERROR("While removing a track. Track collection is empty!");

		TrackCollection::iterator trackIt = std::find( trackCollection->begin() , trackCollection->end() , track );

		if( trackIt != trackCollection->end() ) {
			delete track;
			trackCollection->erase( trackIt );
			return S_OK("Track correctly removed");
		}

		return S_ERROR("While removing a track. Track was not registered by the Track Manager!");
	}


	Return TrackManager::ClearAllContent() {

		if( trackCollection == 0 )
			return S_ERROR("While clearing all content in track manager : assertion trackCollection != 0 failed");

		for( unsigned int i=0 ; i<trackCollection->size() ; i++ ) {
			if( trackCollection->at(i) != 0 )
				delete trackCollection->at(i);
		}
		trackCollection->clear();

		return S_OK("Content cleared in track manager");
	}


	bool TrackManager::TrackContainsHit( Track *track , Hit *hit ) {

		return track->Contains( hit );
	}


}  // namespace 

