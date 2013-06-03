  /// \file Shower.cc
/*
 *
 * Shower.cc source template generated by fclass
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


#include "Objects/Shower.hh"

namespace baboon {

	Shower::Shower() {
		hitToWeightsMap = new HitToWeightsMap;
		startingCone = 0;
		thrust = 0;
		coreCollection = new CoreCollection();
		trackCollection = new TrackCollection();
	}


	Shower::~Shower() {

		hitToWeightsMap->clear();
		delete hitToWeightsMap;
		coreCollection->clear();
		delete coreCollection;
		trackCollection->clear();
		delete trackCollection;
	}


	Return Shower::AddHit( Hit *hit , const double &weight ) {

		if( hit == 0 )
			return BABOON_INVALID_PARAMETER("Assertion hit != 0 failed.");

		if( hitToWeightsMap->empty() ) {
			(*hitToWeightsMap)[ hit ] = weight;
			return BABOON_SUCCESS();
		}

		if( hitToWeightsMap->find( hit ) != hitToWeightsMap->end() )
			return BABOON_ALREADY_PRESENT("Hit is already in the collection");
		else {
			(*hitToWeightsMap)[ hit ] = weight;
			return BABOON_SUCCESS();
		}
	}



	Return Shower::AddHit( Hit *hit ) {

		return AddHit( hit , 1.0 );
	}


	Return Shower::SetHitWeight( Hit *hit , const double &weight ) {

		if( hit == 0 )
			return BABOON_INVALID_PARAMETER("Assertion hit != 0 failed.");

		if( hitToWeightsMap->empty() ) {
			return BABOON_NOT_INITIALIZED("Hit map is empty. Can't set a weight");
		}

		HitToWeightsMap::iterator it = hitToWeightsMap->find( hit );
		if( it != hitToWeightsMap->end() )
			return BABOON_NOT_FOUND("Hit is not in shower. Can't set a weight");
		else {
			it->second = weight;
			return BABOON_SUCCESS();
		}
	}

	Return Shower::RemoveHit( Hit *hit ) {

		if( hitToWeightsMap->empty() )
			return BABOON_NOT_FOUND("Try to remove a hit from an empty hit collection.");

		HitToWeightsMap::iterator it = hitToWeightsMap->find( hit );

		if( it != hitToWeightsMap->end() ) {
			hitToWeightsMap->erase( it );
			return BABOON_SUCCESS();
		}
		else {
			return BABOON_NOT_FOUND("Hit was not in the collection.");
		}

	}

	bool Shower::Contains( Hit *hit ) {

		return ( hitToWeightsMap->find( hit ) != hitToWeightsMap->end() );
	}


	Return Shower::SetStartingPoint( const ThreeVector &startingVec ) {

		startingPoint = startingVec;
		return BABOON_SUCCESS();
	}


	Return Shower::SetStartingCone( Cone *cone ) {

		if( cone == 0 )
			return BABOON_INVALID_PARAMETER("While setting the cone in Shower : assertion cone != 0 failed.");

		if( startingCone != 0 )
			delete startingCone;

		startingCone = cone;
		return BABOON_SUCCESS();
	}


	Return Shower::SetThrust( Track *th ) {

		if( th == 0 )
			return BABOON_INVALID_PARAMETER("Assertion thrust != 0 failed.");
		thrust = th;
		return BABOON_SUCCESS();
	}

	bool Shower::HasThrust() {

		if( thrust == 0 ) return false;
		return true;
	}


	Return Shower::AddCore( Core *core ) {

		if( core == 0 )
			return BABOON_INVALID_PARAMETER("Assertion core != 0 failed.");

		if( std::find( coreCollection->begin() , coreCollection->end() , core ) != coreCollection->end() )
			return BABOON_ALREADY_PRESENT("Core already in the core collection.");

		coreCollection->push_back( core );

		HitCollection *coreHits = core->GetHitCollection();
		for( unsigned int h=0 ; h<coreHits->size() ; h++ ) {
			if( !this->Contains( coreHits->at(h) ) )
				BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , this->AddHit( coreHits->at(h) ) );
		}
		return BABOON_SUCCESS();
	}

	Return Shower::RemoveCore( Core *core ) {

		if( core == 0 )
			return BABOON_INVALID_PARAMETER("Assertion core != 0 failed.");

		CoreCollection::iterator it = std::find( coreCollection->begin() , coreCollection->end() , core );

		if( it != coreCollection->end() ) {
			coreCollection->erase( it );

			HitCollection *coreHits = core->GetHitCollection();
			for( unsigned int h=0 ; h<coreHits->size() ; h++ ) {
				if( this->Contains( coreHits->at(h) ) )
					BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , this->RemoveHit( coreHits->at(h) ) );
			}

			return BABOON_SUCCESS();
		}

		return BABOON_ERROR("Core doesn't exists in core collection.");
	}

	bool Shower::Contains( Core *core ) {

		if( core == 0 ) return false;
		return std::find( coreCollection->begin() , coreCollection->end() , core ) != coreCollection->end();
	}

	Return Shower::AddTrack( Track *track ) {

		if( track == 0 )
			return BABOON_INVALID_PARAMETER("Assertion track != 0 failed.");

		if( std::find( trackCollection->begin() , trackCollection->end() , track ) != trackCollection->end() )
			return BABOON_ALREADY_PRESENT("Track already in the track collection.");

		trackCollection->push_back( track );

		HitCollection *trackHits = track->GetHitCollection();
		for( unsigned int t=0 ; t<trackHits->size() ; t++ ) {
			if( !this->Contains( trackHits->at(t) ) )
				BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , this->AddHit( trackHits->at(t) ) );
		}

		return BABOON_SUCCESS();
	}

	Return Shower::RemoveTrack( Track *track ) {

		if( track == 0 )
			return BABOON_INVALID_PARAMETER("Assertion track != 0 failed.");

		TrackCollection::iterator it = std::find( trackCollection->begin() , trackCollection->end() , track );

		if( it != trackCollection->end() ) {
			trackCollection->erase( it );

			HitCollection *trackHits = track->GetHitCollection();
			for( unsigned int t=0 ; t<trackHits->size() ; t++ ) {
				if( this->Contains( trackHits->at(t) ) )
					BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , this->RemoveHit( trackHits->at(t) ) );
			}

			return BABOON_SUCCESS();
		}

		return BABOON_ERROR("Track doesn't exists in track collection.");
	}

	bool Shower::Contains( Track *track ) {

		if( track == 0 ) return false;
		return std::find( trackCollection->begin() , trackCollection->end() , track ) != trackCollection->end();
	}


}  // namespace 
