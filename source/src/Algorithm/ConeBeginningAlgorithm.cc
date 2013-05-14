  /// \file ConeBeginningAlgorithm.cc
/*
 *
 * ConeBeginningAlgorithm.cc source template generated by fclass
 * Creation date : ven. avr. 19 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Algorithm/ConeBeginningAlgorithm.hh"

using namespace std;

namespace baboon {

	ConeBeginningAlgorithm::ConeBeginningAlgorithm()
		: AbstractAlgorithm("ConeBeginningAlgorithm") {
		needData = false;
	}

	ConeBeginningAlgorithm::~ConeBeginningAlgorithm() {

	}



	void ConeBeginningAlgorithm::Init() {

		hitCollection = HitManager::GetInstance()->GetHitCollection();

	}


	Return ConeBeginningAlgorithm::CheckConsistency() {

		if( hitCollection == 0 )
			return S_ERROR("ConeBeginningAlgorithm bad init. Please check your inputs!");
		return S_OK();
	}


	void ConeBeginningAlgorithm::Execute() {

		vector<IntVector> padsPositionsTemp;
		HitManager* hitMan = HitManager::GetInstance();

		for( unsigned int hitID=0 ; hitID<hitCollection->size() ; hitID++ ) {

			Hit* hit = hitCollection->at(hitID);
			IntVector ijk = hit->GetIJK();


		}

	}


	void ConeBeginningAlgorithm::End() {

	}





}  // namespace 

