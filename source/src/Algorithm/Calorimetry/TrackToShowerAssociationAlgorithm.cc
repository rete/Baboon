  /// \file TrackToShowerAssociationAlgorithm.cc
/*
 *
 * TrackToShowerAssociationAlgorithm.cc source template generated by fclass
 * Creation date : lun. mai 27 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : remi
 */


#include "Algorithm/Calorimetry/TrackToShowerAssociationAlgorithm.hh"

using namespace std;

namespace baboon {

	TrackToShowerAssociationAlgorithm::TrackToShowerAssociationAlgorithm()
		: AbstractAlgorithm("TrackToShowerAssociationAlgorithm") {

		needData = true;
	}

	TrackToShowerAssociationAlgorithm::~TrackToShowerAssociationAlgorithm() {

	}


	Return TrackToShowerAssociationAlgorithm::Init() {

		data.GetValue("Chi2LimitForTrack",&Chi2LimitForTrack);
		data.GetValue("stepLength",&stepLength);
//		data.GetValue("cylinderRadius",&cylinderRadius);
//		data.GetValue("cylinderLength",&cylinderLength);

		return BABOON_SUCCESS();
	}


	Return TrackToShowerAssociationAlgorithm::CheckConsistency() {

		if( ShowerManager::GetInstance()->GetShowerCollection()->empty() )
			return BABOON_NOT_INITIALIZED("No shower stored in the Shower Manager. Please register something before running this!");

		if( TrackManager::GetInstance()->GetTrackCollection()->empty() )
			return BABOON_NOT_INITIALIZED("No track stored in the Track Manager. Please register something before running this!");

		return BABOON_SUCCESS();
	}


	Return TrackToShowerAssociationAlgorithm::Execute() {

/*
		TrackCollection *trackCollection = TrackManager::GetInstance()->GetTrackCollection();
		ShowerCollection *showerCollection = ShowerManager::GetInstance()->GetShowerCollection();
//		HitManager *hitManager = HitManager::GetInstance();

		for( unsigned int i=0 ; i<trackCollection->size() ; i++ ) {

			Track *track = trackCollection->at(i);
			track->SortHits();

			vector<ThreeVector> positions = track->GetIJKs();
			vector<ThreeVector> weigths( track->Size() , ThreeVector(1,1,1) );
			Linear3DFit *fitter = new Linear3DFit(positions,weigths);
			fitter->Fit();

			TrackExtremities extremities = track->GetExtremities();
			ThreeVector firstPoint( extremities.first->GetIJK().at(0) , extremities.first->GetIJK().at(1) , extremities.first->GetIJK().at(2) );
			ThreeVector secondPoint( extremities.second->GetIJK().at(0) , extremities.second->GetIJK().at(1) , extremities.second->GetIJK().at(2) );
			ThreeVector trackDirection = fitter->VectorFromRealLine( firstPoint ) - fitter->VectorFromRealLine( secondPoint );

			if( trackDirection.z() < 0 )
				trackDirection = - trackDirection;


			double normalizedChi2 = fitter->GetChi2()/track->Size();

			cout << "chi2 (normalized) : " << normalizedChi2 << endl;
			cout << "trackDirection : " << trackDirection << endl;

			ThreeVector translation = -trackDirection.unit();

			ThreeVector startingPoint;

			if( firstPoint.z() < secondPoint.z() )
				startingPoint = firstPoint;
			else
				startingPoint = secondPoint;

			bool found = false;
			Shower *showerFound = 0;

			ThreeVector pointedDirection = startingPoint;

			while( true ) {

				pointedDirection += stepLength*translation;

				int I = round( pointedDirection.x() );
				int J = round( pointedDirection.y() );
				int K = round( pointedDirection.z() );

				//  if the point doesn't match to any hit in the calo, it is outside.
				if( !hitManager->PadExists( I , J , K ) )
					break;

				if( !hitManager->PadIsTouched( I , J , K ) )
					continue;

				Hit *hit = hitManager->GetHitAt( I , J , K );

				for( unsigned int s=0 ; s<trackCollection->size() ; s++ ) {

					if( showerCollection->at(s)->Contains( hit ) ) {
						found = true;
						showerFound = showerCollection->at(s);

						cout << "Track association made for shower " << s << endl;
						break;
					}
				}

				if( found ) break;
			}


//			translation.setMag( cylinderLength );
//			ThreeVector x1;
//			ThreeVector x2;

//			if( firstPoint.z() < secondPoint.z() )
//				x1 = firstPoint;
//			else
//				x1 = secondPoint;
//			x2 = x1 - translation;

			if( normalizedChi2 < Chi2LimitForTrack ) {

				Cylinder *cylinder = new Cylinder( x1 , x2 , cylinderRadius );
				int nbOfPads;
				int nbOfFiredPads;

				Shower *predominantShower = 0;
				BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , this->PredominantShowerInCylinder( cylinder , predominantShower ) );

				if( predominantShower == 0 )
					continue;

				for( unsigned int s=0 ; s<showerCollection->size() ; s++ ) {

					Shower *shower = showerCollection->at(s);

					if( shower == predominantShower) {
						shower->AddTrack( track );
						continue;
					}

					else if( shower->Contains( track ) ) {
						BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , shower->RemoveTrack( track ) );
					}
				}
			}

			delete fitter;
		}
		*/
		return BABOON_SUCCESS();
	}


	Return TrackToShowerAssociationAlgorithm::End() {

		return BABOON_SUCCESS();
	}

	Return TrackToShowerAssociationAlgorithm::PredominantShowerInCylinder( Cylinder *cyl , Shower *predominantShower ) {

		/*
		if( cyl == 0 )
			return BABOON_NOT_INITIALIZED("Assertion cylinder != 0 failed");

		HitManager *hitManager = HitManager::GetInstance();

		ThreeVector cylinderCenter = cyl->GetCenter();
		double radius = cyl->GetRadius();
		double length = cyl->GetLength();
		double maximumScale = sqrt( radius*radius + length*length/4.0 );

		HitCollection *hitsInCylinder = new HitCollection();

		for( unsigned int i=-(maximumScale+1) ; i<(maximumScale+1) ; i++ ) {
			for( unsigned int j=-(maximumScale+1) ; j<(maximumScale+1) ; j++ ) {
				for( unsigned int k=-(maximumScale+1) ; k<(maximumScale+1) ; k++ ) {

					if( !hitManager->PadExists( int(cylinderCenter.x())+i , int(cylinderCenter.y())+j , int(cylinderCenter.z())+k ))
						continue;

					if( !hitManager->PadIsTouched( int(cylinderCenter.x())+i , int(cylinderCenter.y())+j , int(cylinderCenter.z())+k ))
						continue;

					Hit *hit = hitManager->GetHitAt( int(cylinderCenter.x())+i , int(cylinderCenter.y())+j , int(cylinderCenter.z())+k );

					if( hit == 0 )
						return BABOON_ERROR("Assertion hit != 0 failed!");

					hitsInCylinder->push_back( hit );
				}
			}
		}

		ShowerManager *showerManager = ShowerManager::GetInstance();
		ShowerCollection *showerCollection = showerManager->GetShowerCollection();

//		IntVector showerCounting(showerCollection->size(),0);
		map<Shower*,int> showerCounting;
//
		for( unsigned int i=0 ; i<showerCollection->size() ; i++ )
			showerCounting[ showerCollection->at(i) ] = 0;
//
//
		for( unsigned int j=0 ; j<hitsInCylinder->size() ; j++ ) {

			Hit *hit = hitsInCylinder->at(j);

			for( unsigned int i=0 ; i<showerCollection->size() ; i++ ) {

				Shower *shower = showerCollection->at(i);

				if( shower->Contains( hit ) )
					showerCounting[ shower ]++;
			}
		}

		map<Shower*,int>::iterator it;
		predominantShower = showerCounting.begin()->first;
		for( it=showerCounting.begin() ; it!=showerCounting.end() ; it++ ) {

			if( it->second > showerCounting[ predominantShower ] )
				predominantShower = it->first;
		}

		cout << "predominantShower : " << predominantShower << endl;


		showerCounting.clear();

		delete hitsInCylinder;

		*/
		return BABOON_SUCCESS();
	}


	Return TrackToShowerAssociationAlgorithm::MergeTrackInShower( Shower *showerToEnlarge , Track *trackToAssociate ) {

		/*
		if( showerToEnlarge == 0 )
			return BABOON_NOT_INITIALIZED("Assertion showerToEnlarge != 0 failed");
		if( trackToAssociate == 0 )
			return BABOON_NOT_INITIALIZED("Assertion trackToAssociate != 0 failed");

		if( showerToEnlarge->IsConnectedTo( trackToAssociate ) )
			return BABOON_SUCCESS();

		TrackCollection *trackCollection = TrackManager::GetInstance()->GetTrackCollection();
		ShowerCollection *showerCollection = ShowerManager::GetInstance()->GetShowerCollection();

		for( unsigned int s=0 ; s<showerCollection->size() ; s++ ) {

			Shower *shower = showerCollection->at(s);




//			for( unsigned int i=0 ; i<trackCollection->size() ; i++ ) {


		}

*/
		return BABOON_SUCCESS();
	}


}  // namespace 
