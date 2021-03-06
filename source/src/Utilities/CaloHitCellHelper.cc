  /// \file HitCellHelper.cc
/*
 *
 * HitCellHelper.cc source template generated by fclass
 * Creation date : ven. juil. 12 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Utilities/CaloHitCellHelper.hh"

using namespace std;

namespace baboon {

	CaloHitCellHelper::CaloHitCellHelper()
		: jFactor(1000),
		  kFactor(1000) {
	}

	CaloHitCellHelper::~CaloHitCellHelper() {

		cellIDToHitMap.clear();
	}


	CaloHit *CaloHitCellHelper::GetCaloHitAt( unsigned int I , unsigned int J , unsigned int K ) {

		unsigned long int key = this->IJKToKey( I , J , K );
		if( cellIDToHitMap.find( key ) == cellIDToHitMap.end() )
			return 0;
		else return cellIDToHitMap[ key ];
	}


	bool CaloHitCellHelper::IsPadFired( unsigned int I , unsigned int J , unsigned int K ) {

		return ( this->GetCaloHitAt( I , J , K ) != 0 );
	}


	unsigned long int CaloHitCellHelper::IJKToKey( unsigned int I , unsigned int J , unsigned int K ) {

		return (kFactor*jFactor)*I + (kFactor)*J + K;
	}

	IntVector CaloHitCellHelper::KeyToIJK( unsigned long int key ) {

		IntVector ijk(3,0);

		ijk.at(0) = ( key / (kFactor*jFactor) );
		ijk.at(1) = ( key % (kFactor*jFactor) / kFactor );
		ijk.at(2) = ( key % kFactor );

		return ijk;
	}



	Return CaloHitCellHelper::AddCaloHit( CaloHit *caloHit ) {

		BABOON_CHECK_POINTER( caloHit );

		IntVector ijk = caloHit->GetIJK();
		unsigned long int key = IJKToKey( ijk.at(0) , ijk.at(1) , ijk.at(2) );

		if( cellIDToHitMap.find( key ) == cellIDToHitMap.end() ) {

			cellIDToHitMap[ key ] = caloHit;
			return BABOON_SUCCESS();
		}
		else
			return BABOON_ALREADY_PRESENT("Hit was already present");
	}


	Return CaloHitCellHelper::RemoveCaloHit( CaloHit *caloHit ) {

		BABOON_CHECK_POINTER( caloHit );

		IntVector ijk = caloHit->GetIJK();
		unsigned long int key = IJKToKey( ijk.at(0) , ijk.at(1) , ijk.at(2) );

		if( cellIDToHitMap.find( key ) == cellIDToHitMap.end() ) {
			return BABOON_NOT_FOUND("Hit not found!");
		}
		else {
			cellIDToHitMap[ key ] = 0;
			cellIDToHitMap.erase( key );
			return BABOON_SUCCESS();
		}
	}



	void CaloHitCellHelper::Clear() {

		cellIDToHitMap.clear();
	}



}  // namespace 

