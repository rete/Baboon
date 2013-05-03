  /// \file NearbyClustering2DAlgorithm.cc
/*
 *
 * NearbyClustering2DAlgorithm.cc source template generated by fclass
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


#include "Algorithm/Clustering/NearbyClustering2DAlgorithm.hh"


using namespace std;

namespace sdhcal {

	NearbyClustering2DAlgorithm::NearbyClustering2DAlgorithm()
		: AbstractAlgorithm("NearbyClustering2DAlgorithm") {
		needData = false;
	}

	NearbyClustering2DAlgorithm::~NearbyClustering2DAlgorithm() {

	}

	void NearbyClustering2DAlgorithm::Init() {

		hitCollection = HitManager::GetInstance()->GetHitCollection();
	}


	Return NearbyClustering2DAlgorithm::CheckConsistency() {

		if( hitCollection == 0 )
			return S_ERROR("NearbyClustering2DAlgorithm bad init. Please check your inputs!");
		return S_OK();
	}


	void NearbyClustering2DAlgorithm::Execute() {

		HitManager *hitManager = HitManager::GetInstance();

		HitCollection tempStorageCollection;

		for( unsigned int hitID=0 ; hitID<hitCollection->size() ; hitID++ ) {

			if( std::find( tempStorageCollection.begin() , tempStorageCollection.end() , hitCollection->at(hitID) ) != tempStorageCollection.end() )
				continue;

			tempStorageCollection.push_back( hitCollection->at(hitID) );

			IntVec ijk1 = hitCollection->at(hitID)->GetIJK();

			for( int i=-1 ; i<=1 ; i++ ) {
				for( int j=-1 ; j<=1 ; j++ ) {

					if( !hitManager->PadExists( ijk1.at(0)+i , ijk1.at(1)+j , ijk1.at(2) ) ) continue;
					if( !hitManager->PadIsTouched( ijk1.at(0)+i , ijk1.at(1)+j , ijk1.at(2) ) ) continue;

					Hit *hit2 = hitManager->GetHitAt( ijk1.at(0)+i , ijk1.at(1)+j , ijk1.at(2) );

					hitCollection->at(hitID)->MergeClusters2D( hit2 );
					tempStorageCollection.push_back( hit2 );

				}
			}




		}

//
//		for( unsigned int hitID=0 ; hitID<hitCollection->size() ; hitID++ ) {
//
//			cout << "cluster size : " << hitCollection->at(hitID)->GetCluster()->GetClusterSize()  << endl;
//
//		}

	}


	void NearbyClustering2DAlgorithm::End() {

	}


}  // namespace 

