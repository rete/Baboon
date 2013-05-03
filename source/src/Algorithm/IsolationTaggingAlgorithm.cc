  /// \file IsolationTaggingAlgorithm.cc
/*
 *
 * IsolationTaggingAlgorithm.cc source template generated by fclass
 * Creation date : lun. avr. 22 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */




#include "Algorithm/IsolationTaggingAlgorithm.hh"



using namespace std;



namespace baboon {


	IsolationTaggingAlgorithm::IsolationTaggingAlgorithm()
		: AbstractAlgorithm("IsolationTaggingAlgorithm") {
		needData = false;
		hitCollection = 0;
		// default values
		distance = 2;
		concentrationLimit = 0.2;
	}


	IsolationTaggingAlgorithm::~IsolationTaggingAlgorithm() {

	}


	void IsolationTaggingAlgorithm::Init() {

		data.GetValue("distance",&distance);
		data.GetValue("concentrationLimit",&concentrationLimit);
		hitCollection = HitManager::GetInstance()->GetHitCollection();
		isolationWeights.clear();

	}

	Return IsolationTaggingAlgorithm::CheckConsistency() {

		if( hitCollection == 0 )
			return S_ERROR("IsolationTaggingAlgorithm bad init. Please check your inputs!");
		return S_OK();
	}

	void IsolationTaggingAlgorithm::Execute() {

		HitManager *hitManager = HitManager::GetInstance();

		unsigned int hitID = 0;
		unsigned int size = hitCollection->size();
		unsigned int volume = (2*distance+1)*(2*distance+1)*(2*distance+1);

		for( hitID=0 ; hitID<size ; hitID++ ) {

			Hit *hit = hitCollection->at(hitID);
			IntVec ijk1 = hit->GetIJK();

			int hitCount = 0;

			for( int i=-distance ; i<=distance ; i++ ) {
				for( int j=-distance ; j<=distance ; j++ ) {
					for( int k=-distance ; k<=distance ; k++ ) {

						if( !hitManager->PadExists( ijk1.at(0)+i , ijk1.at(1)+j , ijk1.at(2)+k ) ) continue;
						if( !hitManager->PadIsTouched( ijk1.at(0)+i , ijk1.at(1)+j , ijk1.at(2)+k ) ) continue;

						hitCount++;
//						Hit *hit2 = hitManager->GetHitAt( ijk1.at(0)+i , ijk1.at(1)+j , ijk1.at(2)+k );

					}
				}
			}
			isolationWeights.push_back( double(hitCount) / volume );

			int clusterSize = hit->GetCluster3D()->GetClusterSize();

//			if( clusterSize < 4 ) {
			if( double(hitCount) / volume < concentrationLimit ) {
				hit->SetHitTag( fIsolated );
			}
		}
	}



	void IsolationTaggingAlgorithm::End() {
		hitCollection = 0;
	}





}  // namespace 

