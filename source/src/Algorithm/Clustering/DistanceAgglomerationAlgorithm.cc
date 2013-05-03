  /// \file DistanceAgglomerationAlgorithm.cc
/*
 *
 * DistanceAgglomerationAlgorithm.cc source template generated by fclass
 * Creation date : mer. avr. 17 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */



#include "Algorithm/Clustering/DistanceAgglomerationAlgorithm.hh"


namespace sdhcal {

	DistanceAgglomerationAlgorithm::DistanceAgglomerationAlgorithm()
		: AbstractAlgorithm("DistanceAgglomerationAlgorithm") {
		hitCollection = 0;
		startingHit = 0;
		needData = false;
	}

	DistanceAgglomerationAlgorithm::~DistanceAgglomerationAlgorithm() {

	}


	void DistanceAgglomerationAlgorithm::Init() {

		if( aggregate != 0 ) delete aggregate;
		aggregate = new Cluster();

	}


	Return DistanceAgglomerationAlgorithm::CheckConsistency() {

		if (hitCollection == 0 || startingHit == 0)
			return S_ERROR("DistanceAgglomerationAlgorithm bad init. Please check your inputs!");
		if( aggregate != 0 )
			if( aggregate->GetClusterSize() != 0 )
				return S_ERROR("DistanceAgglomerationAlgorithm bad check on cluster size (!= 0). Bad init of aggregate pointer!");

		return S_OK();
	}


	void DistanceAgglomerationAlgorithm::Execute() {

//		for( unsigned int i=0 ; i<hitCollection->size() ; i++ ) {
//
//
//
//		}

	}


	void DistanceAgglomerationAlgorithm::End() {

	}


}  // namespace 
