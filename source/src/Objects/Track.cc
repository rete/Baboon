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



#include "Objects/Track.hh"


using namespace std;


namespace baboon {

	Track::Track()
	: HitCompositeObject(),
	  TypedObject("Track") {

	}

	Track::~Track() {

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


	std::vector<ThreeVector> Track::GetPositions() const {

		std::vector<ThreeVector> vecCol;
		for( unsigned int i=0 ; i<hitCollection->size() ; i++ ) {
			vecCol.push_back( hitCollection->at(i)->GetPosition() );
		}
		return vecCol;
	}

	std::vector<ThreeVector> Track::GetIJKs() const {

		std::vector<ThreeVector> ijks;
		for( unsigned int i=0 ; i<hitCollection->size() ; i++ ) {
			ThreeVector ijk( hitCollection->at(i)->GetIJK().at(0) , hitCollection->at(i)->GetIJK().at(1) , hitCollection->at(i)->GetIJK().at(2) );
			ijks.push_back( ijk );
		}
		return ijks;
	}

	Return Track::SortHits() {

		if( hitCollection->size() <= 1 )
			return BABOON_SUCCESS("Sort one element or an empty list is easy...");

		int i = 0;
		int j = 0;
		Hit *hit = 0;

		for( j=1 ; j<hitCollection->size() ; j++ ) {

			i = j-1;
			while( hitCollection->at(j)->GetIJK().at(2) < hitCollection->at(i)->GetIJK().at(2) ) {
				hit = hitCollection->at(i);
				hitCollection->at(i) = hitCollection->at(j);
				hitCollection->at(j) = hit;
				i=i-1;
				j=j-1;
				if( i<0 ) break;
			}
		}
		return BABOON_SUCCESS();
	}




}  // namespace 

