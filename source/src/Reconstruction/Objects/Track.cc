  /// \file Track.cc
/*
 *
 * Track.cc source template generated by fclass
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



#include "Reconstruction/Objects/Track.hh"


using namespace std;


namespace baboon {

	Track::Track()
		: RecoObject() {
		hitCollection = new HitCollection;
	}

	Track::~Track() {
		delete hitCollection;
	}

	TrackExtremities Track::GetExtremities() {

		TrackExtremities extrem;

		bool firstIt = true;
		Hit *hitFirst = 0;
		Hit *hitSecond = 0;
		for(unsigned int i=0 ; i<hitCollection->size() ; i++) {
			Hit *hit1 = hitCollection->at(i);
			ThreeVector vec1 = hit1->GetPosition();
			for(unsigned int j=0 ; j<hitCollection->size() ; j++) {
				if(i == j) continue;
				Hit *hit2 = hitCollection->at(j);
				ThreeVector vec2 = hit2->GetPosition();
				if(firstIt) {
					hitFirst = hit1;
					hitSecond = hit2;
//					extrem = make_pair( hit1 , hit2 );
					firstIt = false;
					continue;
				}
				// compare the distance between the two hits with the current minimum distance (in extrem)
				if( (hit1->GetPosition() - hit2->GetPosition()).mag()
						> (hitFirst->GetPosition() - hitSecond->GetPosition()).mag() ) {
					hitFirst = hit1;
					hitSecond = hit2;
				}
			}
		}

		extrem = make_pair( hitFirst , hitSecond );

		return extrem;
	}


	double Track::GetTrackLength() {
		ThreeVector vec( extremities.first->GetPosition() - extremities.second->GetPosition() );
		return vec.mag();
	}


	std::vector<ThreeVector> Track::GetPositions() const {

		std::vector<ThreeVector> vecCol;
		for( unsigned int i=0 ; i<hitCollection->size() ; i++ ) {
			vecCol.push_back( hitCollection->at(i)->GetPosition() );
		}
		return vecCol;

	}

	Return Track::AddHit( Hit *hit ) {

		if( find( hitCollection->begin() , hitCollection->end() , hit ) == hitCollection->end() )
			hitCollection->push_back(hit);
		return S_OK();
	}

	Return Track::RemoveHit( Hit *hit ) {

		HitCollection::iterator pos = find( hitCollection->begin() , hitCollection->end() , hit );
		if( pos == hitCollection->end() )
			return S_ERROR("Hit was not in track");
		else {
			hitCollection->erase(pos);
			return S_OK("Hit removed from a track");
		}
	}


	bool Track::Contains( Hit *hit ) {
		return ( std::find( hitCollection->begin() , hitCollection->end() , hit ) != hitCollection->end() );
	}


	Return Track::OrderHits() {



	}




}  // namespace 

