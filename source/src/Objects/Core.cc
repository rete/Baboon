  /// \file Core.cc
/*
 *
 * Core.cc source template generated by fclass
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


#include "Objects/Core.hh"

using namespace std;

namespace baboon {


	Core::Core() {

		hitCollection = new HitCollection();
	}

	Core::~Core() {
		delete hitCollection;

	}

	Return Core::AddHit( Hit *hit ) {

		if( hit == 0 )
			return BABOON_INVALID_PARAMETER("Assertion hit != 0 failed");

		if( hitCollection == 0 )
			return BABOON_NOT_INITIALIZED("Hit collection is a null pointer. Something very bad happened...");

		if( hitCollection->empty() ) {
			hitCollection->push_back( hit );
			return BABOON_SUCCESS();
		}

		if( std::find( hitCollection->begin() , hitCollection->end() , hit ) != hitCollection->end() )
			return BABOON_ALREADY_PRESENT("Hit was already in the collection");
		else {
			hitCollection->push_back( hit );
			return BABOON_SUCCESS();
		}
	}


	Return Core::RemoveHit( Hit *hit ) {

		if( hitCollection == 0 )
			return BABOON_INVALID_PARAMETER("Hit collection is a null pointer. Something very bad happened...");

		if( hitCollection->empty() )
			return BABOON_ALREADY_PRESENT("Try to remove a hit from an empty hit collection.");

		HitCollection::iterator it = std::find( hitCollection->begin() , hitCollection->end() , hit );

		if( it != hitCollection->end() ) {
			hitCollection->erase( it );
			return BABOON_SUCCESS();
		}
		else {
			return BABOON_NOT_FOUND("Hit was not in the collection.");
		}
	}

	bool Core::Contains( Hit *hit ) {

		return ( find( hitCollection->begin() , hitCollection->end() , hit ) != hitCollection->end() );
	}


	Return Core::SetStartingCone( Cone *cone ) {

		if( cone == 0 )
			return BABOON_INVALID_PARAMETER("Assertion cone != 0 failed");
		startingCone = cone;
		return BABOON_SUCCESS();
	}


	ThreeVector Core::CenterOfGravity() {

		return hitCollection->GetBarycenter();
	}


	int Core::GetFirstHitLayer() {

		int firstLayer = 10000;
		for( unsigned int i=0 ; i<hitCollection->size() ; i++ ) {
			if( hitCollection->at(i)->GetIJK().at(2) < firstLayer ) {
				firstLayer = hitCollection->at(i)->GetIJK().at(2);
			}
		}

		return firstLayer;
	}


	int Core::GetLastHitLayer() {

		int lastLayer = 0;
		for( unsigned int i=0 ; i<hitCollection->size() ; i++ )
			if( hitCollection->at(i)->GetIJK().at(2) > lastLayer )
				lastLayer = hitCollection->at(i)->GetIJK().at(2);

		return lastLayer;
	}

}  // namespace 
