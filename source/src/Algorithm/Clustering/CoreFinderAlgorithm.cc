/*
 *
 * CoreFinderAlgorithm.cc cpp file template generated by fclass
 * Creation date : Mon Apr  8 17:46:03 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#include "Algorithm/Clustering/CoreFinderAlgorithm.hh"


using namespace std ;

namespace baboon {


	CoreFinderAlgorithm::CoreFinderAlgorithm()
		: AbstractAlgorithm("CoreFinderAlgorithm") {
		hitCollection = 0;
		distance = 1;
		needData = true;
	}

	CoreFinderAlgorithm::~CoreFinderAlgorithm() {}


	Return CoreFinderAlgorithm::Init() {

		hitCollection = HitManager::GetInstance()->GetHitCollection();
		threshCountVec.clear();
		data.GetValue("coreCountThreshold",&coreCountThreshold);
		data.GetValue("minimumThresholdConcentration",&minimumThresholdConcentration);
		return BABOON_SUCCESS();
	}

	Return CoreFinderAlgorithm::CheckConsistency() {

		if(hitCollection == 0)
			return BABOON_ERROR("CoreFinderAlgorithm bad init. Please check your inputs!");
		return BABOON_SUCCESS();

	}

	Return CoreFinderAlgorithm::Execute() {

		HitManager *hitManager = HitManager::GetInstance();
		int nbOfCoreHits = 0;
		for( unsigned int l=0 ; l<hitCollection->size() ; l++ ) {

			Hit *hit = hitCollection->at(l);
			IntVector ijk = hit->GetIJK();
			HitThreshold thresh = hit->GetThreshold();


			int count = 0;

			for(int i=-1 ; i<=1 ; i++ ) {
				for(int j=-1 ; j<=1 ; j++) {
					for(int k=-1 ; k<=1 ; k++) {

						if( !hitManager->PadExists( ijk.at(0)+i , ijk.at(1)+j , ijk.at(2)+k ) ) continue;
						if( hitManager->PadIsTouched( ijk.at(0)+i , ijk.at(1)+j , ijk.at(2)+k ) ) {
							Hit* touchedHit = hitManager->GetHitAt( ijk.at(0)+i , ijk.at(1)+j , ijk.at(2)+k );
							int factor = 1;
							factor *= (1+touchedHit->GetThreshold());
							count += factor;
						}

					}
				}
			}

			if(count / 27.0 >= minimumThresholdConcentration ) {
				hit->SetHitTag( fCore );
				nbOfCoreHits++;
			}
		}
		return BABOON_SUCCESS();
	}

	Return CoreFinderAlgorithm::End() {

		return BABOON_SUCCESS();
	}

}
