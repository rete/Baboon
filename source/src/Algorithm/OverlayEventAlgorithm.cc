  /// \file OverlayEventAlgorithm.cc
/*
 *
 * OverlayEventAlgorithm.cc source template generated by fclass
 * Creation date : jeu. nov. 7 2013
 *
 * This file is part of XXX libraries.
 * 
 * XXX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * XXX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author : R�mi Et�
 * @version
 * @copyright
 *
 *
 */


#include "Algorithm/OverlayEventAlgorithm.hh"


#include "Utilities/CaloHitHelper.hh"
#include "Managers/AlgorithmManager.hh"
#include "Algorithm/Calorimetry/TrackFinderAlgorithm.hh"
#include "Algorithm/Calorimetry/ClusteringAlgorithm.hh"
#include "Managers/TrackManager.hh"
#include "Managers/ClusteringManager.hh"
#include "Managers/ExternInfoManager.hh"

namespace baboon {

	OverlayEventAlgorithm::OverlayEventAlgorithm()
		: AbstractAlgorithm("OverlayEventAlgorithm") {

		needData = true;
		_overlayDone = false;
		_generatesLCTracks = false;
		_collection1 = nullptr;
		_collection2 = nullptr;
		_collectionToOverlay1 = new CaloHitCollection();
		_collectionToOverlay2 = new CaloHitCollection();
		_lostHitCollection = new CaloHitCollection();
		_overlaidCollection = new CaloHitCollection();
		_calorimeter1 = nullptr;
		_calorimeter2 = nullptr;
	}

	OverlayEventAlgorithm::~OverlayEventAlgorithm() {

	}


	Return OverlayEventAlgorithm::Init() {

		data.GetValue( "useTrackInfo" , &_useTrackInfo );
		data.GetValue( "generatesLCTracks" , &_generatesLCTracks );
		data.GetValue( "separationDistance" , &_separationDistance );
		data.GetValue( "particleType1" , &_particleType1 );
		data.GetValue( "particleType2" , &_particleType2 );
		data.GetValue( "inputEnergy1" , &_inputEnergy1 );
		data.GetValue( "inputEnergy2" , &_inputEnergy2 );

		_lostHitCollection->clear();
		_collectionToOverlay1->clear();
		_collectionToOverlay2->clear();
		_overlaidCollection->clear();

		_overlayDone = false;

		_trackPair.first = 0;
		_trackPair.second = 0;

		return BABOON_SUCCESS();
	}


	Return OverlayEventAlgorithm::CheckConsistency() {

		BABOON_CHECK_POINTER( _calorimeter1 );
		BABOON_CHECK_POINTER( _calorimeter2 );
		BABOON_CHECK_POINTER( _collection1 );
		BABOON_CHECK_POINTER( _collection2 );

		if( _particleType1 != "charged" && _particleType1 != "neutral" )
			return BABOON_ERROR("OverlayEventAlgortihm : particle type 1 has a bad type");

		if( _particleType2 != "charged" && _particleType2 != "neutral" )
			return BABOON_ERROR("OverlayEventAlgortihm : particle type 2 has a bad type");

		return BABOON_SUCCESS();
	}


	Return OverlayEventAlgorithm::Execute() {

		if( _useTrackInfo ) {

			ExternInfoManager *extInfoMgr = ExternInfoManager::GetInstance();
			AlgorithmManager *algoMan = AlgorithmManager::GetInstance();
			TrackFinderAlgorithm *trackFinder( 0 );
			ClusteringAlgorithm *clusteringAlgo( 0 );

			if( algoMan->AlgorithmIsRegistered("TrackFinderAlgorithm") ) {
				trackFinder = (TrackFinderAlgorithm *) algoMan->GetAlgorithm("TrackFinderAlgorithm");
			}
			else
				return BABOON_ERROR("TrackFinder algo was not registered for OverlayEvent algo. "
						"Please register it before running or disable the \"useTrackInfo\" parameter");

			if( algoMan->AlgorithmIsRegistered("ClusteringAlgorithm") ) {
				clusteringAlgo = (ClusteringAlgorithm *) algoMan->GetAlgorithm("ClusteringAlgorithm");
			}
			else
				return BABOON_ERROR("ClusteringAlgorithm algo was not registered for OverlayEvent algo. "
						"Please register it before running or disable the \"useTrackInfo\" parameter");

			ThreeVector trackEndPosition1;
			ThreeVector trackEndPosition2;
			ThreeVector trackBeginPosition1;
			ThreeVector trackBeginPosition2;
			ThreeVector backwardThrust1;
			ThreeVector backwardThrust2;

//----------------------------------------------------------------------------------------------------

			for( unsigned int h=0 ; h<_collection1->size() ; h++ )
				_calorimeter1->AddCaloHit( _collection1->at(h) );


			ClusterCollection *clusters = new ClusterCollection();
			clusteringAlgo->SetClusteringMode( fClustering2D );
			clusteringAlgo->SetCalorimeter( _calorimeter1 );
			clusteringAlgo->SetClusterCollection( clusters );
			clusteringAlgo->Process();

			// Register them to the clustering manager just for memory management...
			for( unsigned int c=0 ; c<clusters->size() ; c++ ) {
				BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , ClusteringManager::GetInstance()->AddCluster( clusters->at(c) ) );
			}

			clusters->clear();

			// Perform track finder
			trackFinder->SetCalorimeter( _calorimeter1 );
			trackFinder->Process();

			TrackManager *trackMan = TrackManager::GetInstance();
			TrackCollection *trackCollection1 = trackMan->GetTrackCollection();
			ThreeVector showerEnteringPoint1(this->FindShowerEnteringPoint(_collection1));
			bool primaryTrackFound1 = false;

			for( unsigned int tr=0 ; tr<trackCollection1->size() ; tr++ ) {

				Track *track1 = trackCollection1->at(tr);
				OverlayEventAlgorithm::TrackInfo *trackInfo = new OverlayEventAlgorithm::TrackInfo;
				trackInfo->track = track1;
				this->FillTrackInfo( trackInfo );

				if( trackInfo->isPrimaryTrack ) {

					primaryTrackFound1 = true;
					if( _particleType1 == "neutral" )
						this->EraseTrackFromCollection( trackInfo , _collection1 );

					trackEndPosition1 = trackInfo->endPosition;
					trackBeginPosition1 = trackInfo->beginPosition;
					backwardThrust1 = trackInfo->backwardThrust;

					if( _particleType1 == "charged" ) {
						baboon::TrackInfo *externTrackInfo = extInfoMgr->CreateTrackInfo();

						ThreeVector translation(_calorimeter1->GetRepeatX()/2.0 + _separationDistance / 2.0 - showerEnteringPoint1.x()
								, _calorimeter1->GetRepeatY()/2.0 - showerEnteringPoint1.y()
								, 0.0);
						externTrackInfo->enteringPoint = showerEnteringPoint1 + translation;
						externTrackInfo->momentum = backwardThrust1*_inputEnergy1;
						externTrackInfo->charge = 1;
				 }

					delete trackInfo;
					break;
				}

				delete trackInfo;
			}

			trackMan->ClearAllContent();
			ClusteringManager::GetInstance()->ClearAllContent();

			// If primary track is not found and we asked for a charged particle, the overlay is not done.
			if( !primaryTrackFound1 && _particleType1 == "neutral" ) {

				delete clusters;
				clusters = 0;
				return BABOON_SUCCESS();
			}

			ThreeVector translation1;

//			if( _particleType1 == "neutral" ) {
//
//				ThreeVector cog1;
//				BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , CaloHitHelper::ComputeBarycenter( _collection1 , cog1 ) );
//				translation1 = ThreeVector( _calorimeter1->GetRepeatX()/2.0 + _separationDistance / 2.0 - cog1.x()
//										, _calorimeter1->GetRepeatY()/2.0 - cog1.y()
//										, 0.0 );
//			}
//			else {
				translation1 = ThreeVector( _calorimeter1->GetRepeatX()/2.0 + _separationDistance / 2.0 - showerEnteringPoint1.x()
										, _calorimeter1->GetRepeatY()/2.0 - showerEnteringPoint1.y()
										, 0.0 );

//			}

//----------------------------------------------------------------------------------------------------

			for( unsigned int h=0 ; h<_collection2->size() ; h++ )
				_calorimeter2->AddCaloHit( _collection2->at(h) );


			clusteringAlgo->SetClusteringMode( fClustering2D );
			clusteringAlgo->SetCalorimeter( _calorimeter2 );
			clusteringAlgo->SetClusterCollection( clusters );
			clusteringAlgo->Process();

			// Register them to the clustering manager just for memory management...
			for( unsigned int c=0 ; c<clusters->size() ; c++ ) {
				BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , ClusteringManager::GetInstance()->AddCluster( clusters->at(c) ) );
			}

			clusters->clear();

			// Perform track finder
			trackFinder->SetCalorimeter( _calorimeter2 );
			trackFinder->Process();

			TrackCollection *trackCollection2 = trackMan->GetTrackCollection();
			ThreeVector showerEnteringPoint2(this->FindShowerEnteringPoint(_collection2));
			bool primaryTrackFound2 = false;

			for( unsigned int tr=0 ; tr<trackCollection2->size() ; tr++ ) {

				Track *track2 = trackCollection2->at(tr);
				TrackInfo *trackInfo = new TrackInfo;
				trackInfo->track = track2;
				this->FillTrackInfo( trackInfo );

				if( trackInfo->isPrimaryTrack ) {

					primaryTrackFound2 = true;
					if( _particleType2 == "neutral" )
						this->EraseTrackFromCollection( trackInfo , _collection2 );

					trackEndPosition2 = trackInfo->endPosition;
					trackBeginPosition2 = trackInfo->beginPosition;
					backwardThrust2 = trackInfo->backwardThrust;

					if( _particleType2 == "charged" ) {
						baboon::TrackInfo *externTrackInfo = extInfoMgr->CreateTrackInfo();
						ThreeVector translation( _calorimeter2->GetRepeatX()/2.0 - _separationDistance / 2.0 - showerEnteringPoint2.x()
								, _calorimeter2->GetRepeatY()/2.0 - showerEnteringPoint2.y()
								, 0.0 );
						externTrackInfo->enteringPoint = showerEnteringPoint2 + translation;
						externTrackInfo->momentum = backwardThrust2*_inputEnergy2;
						externTrackInfo->charge = 1;
				 }

					delete trackInfo;
					break;
				}

				delete trackInfo;
			}

			trackMan->ClearAllContent();
			ClusteringManager::GetInstance()->ClearAllContent();

			delete clusters;
			clusters = 0;

			if( !primaryTrackFound2  && _particleType2 == "neutral" ) {

				return BABOON_SUCCESS();
			}

			ThreeVector translation2;

//			if( _particleType2 == "neutral" ) {
//
//				ThreeVector cog2;
//				BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , CaloHitHelper::ComputeBarycenter( _collection2 , cog2 ) );
//				translation2 = ThreeVector( _calorimeter2->GetRepeatX()/2.0 - _separationDistance / 2.0 - cog2.x()
//										, _calorimeter2->GetRepeatY()/2.0 - cog2.y()
//										, 0.0 );
//			}
//			else {
				translation2 = ThreeVector( _calorimeter2->GetRepeatX()/2.0 - _separationDistance / 2.0 - showerEnteringPoint2.x()
										, _calorimeter2->GetRepeatY()/2.0 - showerEnteringPoint2.y()
										, 0.0 );

//			}

			this->TranslateCollection( _calorimeter1 , _collection1 , _collectionToOverlay1 , translation1 );
			this->TranslateCollection( _calorimeter2 , _collection2 , _collectionToOverlay2 , translation2 );

			this->OverlayCollections( _collectionToOverlay1 , _collectionToOverlay2 );

			_overlayDone = true;

		}
		else {

			// Overlay the 2 collections wrt the barycenter.

			ThreeVector cog1;
			ThreeVector cog2;

			BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , CaloHitHelper::ComputeBarycenter( _collection1 , cog1 ) );
			BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , CaloHitHelper::ComputeBarycenter( _collection2 , cog2 ) );

			ThreeVector translation1( _calorimeter1->GetRepeatX()/2.0 + _separationDistance / 2.0 - cog1.x()
									, _calorimeter1->GetRepeatY()/2.0 - cog1.y()
									, 0.0 );
			ThreeVector translation2( _calorimeter2->GetRepeatX()/2.0 - _separationDistance / 2.0 - cog2.x()
									, _calorimeter2->GetRepeatY()/2.0 - cog2.y()
									, 0 );

			this->TranslateCollection( _calorimeter1 , _collection1 , _collectionToOverlay1 , translation1 );
			this->TranslateCollection( _calorimeter2 , _collection2 , _collectionToOverlay2 , translation2 );

			this->OverlayCollections( _collectionToOverlay1 , _collectionToOverlay2 );

			_overlayDone = true;
		}



		return BABOON_SUCCESS();
	}


	Return OverlayEventAlgorithm::End() {

		_calorimeter1 = 0;
		_calorimeter2 = 0;

		return BABOON_SUCCESS();
	}




	void OverlayEventAlgorithm::TranslateCollection( Calorimeter *_calorimeter
													  , CaloHitCollection *initialCollection
													  , CaloHitCollection *finalCollection
													  , const ThreeVector &vec ) {

		if( initialCollection == 0 || finalCollection == 0 )
			return;

		for( unsigned int h=0 ; h<initialCollection->size() ; h++ ) {

			CaloHit *caloHit = initialCollection->at(h);
			IntVector ijk = caloHit->GetIJK();

			// if the hit is outside of the sdhcal remove it from the collection!
			if( ( ijk.at(0) + round( vec.x() ) ) < 0
			||  ( ijk.at(1) + round( vec.y() ) ) < 0
			||  ( ijk.at(2) + round( vec.z() ) ) < 0
			||  ( ijk.at(0) + round( vec.x() ) ) > _calorimeter->GetRepeatX()
			||  ( ijk.at(1) + round( vec.y() ) ) > _calorimeter->GetRepeatY()
			||  ( ijk.at(2) + round( vec.z() ) ) > _calorimeter->GetNbOfLayers() ) {

				_lostHitCollection->push_back( caloHit );
			}
			else {
				caloHit->SetIJK( ijk.at(0) + round( vec.x() )
								,ijk.at(1) + round( vec.y() )
								,ijk.at(2) + round( vec.z() ) );

				ThreeVector pos;
				double xShift = _calorimeter->GetCellSize0();
				double yShift = _calorimeter->GetCellSize1();
				double zShift = _calorimeter->GetLayerThickness();
				pos.setX( caloHit->GetPosition().x() + xShift*round( vec.x()) );
				pos.setY( caloHit->GetPosition().y() + yShift*round( vec.y()) );
				pos.setZ( caloHit->GetPosition().z() + zShift*round( vec.z()) );
				caloHit->SetPosition(pos);

				finalCollection->push_back( caloHit );
			}
		}
	}



	void OverlayEventAlgorithm::OverlayCollections( CaloHitCollection *collection1 , CaloHitCollection *collection2 ) {

		if( collection1 == 0 || collection2 == 0 )
			return;

		for( unsigned int h=0 ; h<collection1->size() ; h++ )
			_overlaidCollection->push_back( collection1->at(h) );


		int count = 0;
		CaloHitThreshold thresholdTable[3][3] = { fCaloHitThr1 , fCaloHitThr2 , fCaloHitThr3 ,
												  fCaloHitThr2 , fCaloHitThr2 , fCaloHitThr3 ,
												  fCaloHitThr3 , fCaloHitThr3 , fCaloHitThr3 };

		for( unsigned int h2=0 ; h2<collection2->size() ; h2++ ) {

			CaloHit *caloHit2 = collection2->at(h2);
			IntVector ijk2 = caloHit2->GetIJK();

			bool hitIsOverlaid = false;

			for( unsigned int h1=0 ; h1<collection1->size() ; h1++ ) {

				CaloHit *caloHit1 = collection1->at(h1);
				IntVector ijk1 = caloHit1->GetIJK();

				if( ijk1.at(0) == ijk2.at(0)
				 && ijk1.at(1) == ijk2.at(1)
				 && ijk1.at(2) == ijk2.at(2) ) {

					hitIsOverlaid = true;
					count ++;
					caloHit1->SetTypeID( 3 );    // Overlaid hit type

					CaloHitThreshold fThr1 = caloHit1->GetThreshold();
					CaloHitThreshold fThr2 = caloHit2->GetThreshold();
					CaloHitThreshold newThr = thresholdTable [ThresholdToInt(fThr1)][ThresholdToInt(fThr2)];

					if( newThr == fThr1 ) continue;
					else {
						caloHit1->SetThreshold( newThr );
					}

				}
				if( hitIsOverlaid ) break;
			}

			if( !hitIsOverlaid ) _overlaidCollection->push_back( caloHit2 );

		}

	}

	int OverlayEventAlgorithm::ThresholdToInt( CaloHitThreshold fThr ) {

		if( fThr == fCaloHitThr1 ) return 0;
		else if( fThr == fCaloHitThr2 ) return 1;
		else return 2;
	}


	void OverlayEventAlgorithm::FillTrackInfo( TrackInfo *trackInfo ) {


		bool isPrimaryTrack = false;
		Track *track = trackInfo->track;
		trackInfo->isPrimaryTrack = false;
		track->SortHits();
		CaloHitCollection *trackHits = track->GetCaloHitCollection();
		int lastTrackLayer = track->GetCaloHitCollection()->at( track->Size() - 1 )->GetIJK().at(2);
		int firstTrackLayer = track->GetCaloHitCollection()->at( 0 )->GetIJK().at(2);

		ClusterCollection trackClusters;
		Cluster *currentCluster = new Cluster();

		for( unsigned int h=0 ; h<trackHits->size() ; h++ ) {

			CaloHit *caloHit = trackHits->at( h );
			if( h == trackHits->size() - 1 ) {
				currentCluster->AddCaloHit( caloHit );
				trackClusters.push_back( currentCluster );
				break;
			}

			CaloHit *hit = trackHits->at( h );
			CaloHit *nextHit = trackHits->at( h+1 );

			if( hit->GetIJK().at(2) == nextHit->GetIJK().at(2) ) {
				currentCluster->AddCaloHit( hit );
				continue;
			}
			else {
				currentCluster->AddCaloHit( hit );
				trackClusters.push_back( currentCluster );
				currentCluster = new Cluster();
			}
		}

		currentCluster = 0;
		Cluster *firstCluster( 0 );
		ThreeVector backwardThrust;

		for( unsigned int cl=0 ; cl<trackClusters.size() ; cl++ ) {

			if( cl == 0 ) {
				firstCluster = trackClusters.at( cl );
				trackInfo->beginPosition = trackClusters.at( cl )->GetPosition( fComputeCell );
			}

			if( cl == trackClusters.size() - 1 ) {

				trackInfo->endPosition = trackClusters.at( cl )->GetPosition( fComputeCell );
				break;
			}

			Cluster *cluster = trackClusters.at( cl );
			Cluster *nextCluster = trackClusters.at( cl+1 );
			ThreeVector difference = nextCluster->GetPosition( fComputePosition ) - cluster->GetPosition( fComputePosition );

			if( difference.z() < 0 )
				difference = -difference;

			// for the track beginning
			if( cl < 4 ) {
				backwardThrust += difference;

			}
		}

		if( backwardThrust != ThreeVector() )
			backwardThrust.setMag( 1.f );


		if( backwardThrust.z() < 0 )
			backwardThrust = -backwardThrust;

		trackInfo->backwardThrust = backwardThrust;

		if( backwardThrust.theta() < 0.3 ) {
			if( firstCluster->GetPosition( fComputeCell ).z() < 3 ) {
				trackInfo->isPrimaryTrack = true;
			}
		}

		for( unsigned int cl=0 ; cl<trackClusters.size() ; cl++ ) {
			if( trackClusters.at( cl ) != 0 )
				delete trackClusters.at( cl );
		}
		trackClusters.clear();

	}



	void OverlayEventAlgorithm::EraseTrackFromCollection( Track *track , CaloHitCollection *collection ) {

		if( track == 0 || collection == 0 )
			return;

		CaloHitCollection *trackHits = track->GetCaloHitCollection();

		for( unsigned int h=0 ; h<collection->size() ; h++ ) {

			CaloHit *caloHit = collection->at( h );

			for( unsigned int trH=0 ; trH<trackHits->size() ; trH++ ) {

				CaloHit *trackHit = trackHits->at(trH);
				if( caloHit == trackHit ) {

					collection->erase( collection->begin()  + h );
					h--;
					break;
				}

			}
		}

	}


	void OverlayEventAlgorithm::EraseTrackFromCollection( TrackInfo *trackInfo , CaloHitCollection *collection )
	{
		unsigned int layerCut = round(trackInfo->endPosition.z());

		for(unsigned int i=0 ; i<collection->size() ; i++)
		{
			if(collection->at(i)->GetIJK().at(2) < layerCut)
			{
				collection->erase( collection->begin()  + i );
				i--;
			}
		}
	}


	ThreeVector OverlayEventAlgorithm::FindShowerEnteringPoint(CaloHitCollection *collection)
	{
		float xyCut = 4.0;
		unsigned int layerCut = 4;
		float barycenterX = 0.f;
		float barycenterY = 0.f;

		for(unsigned int c=0 ; c<collection->size() ; c++)
		{
			CaloHit *pCaloHit = collection->at(c);

			barycenterX += pCaloHit->GetIJK().at(0);
			barycenterY += pCaloHit->GetIJK().at(1);
		}

		barycenterX /= collection->size();
		barycenterY /= collection->size();

		CaloHitCollection caloHitVecForEnteringPoint;

		for(unsigned int c=0 ; c<collection->size() ; c++)
		{
			CaloHit *pCaloHit = collection->at(c);

			if(fabs(pCaloHit->GetIJK().at(0)-barycenterX) < xyCut
					&& fabs(pCaloHit->GetIJK().at(1)-barycenterY) < xyCut
					&& pCaloHit->GetIJK().at(2) < layerCut)
			{
				caloHitVecForEnteringPoint.push_back(pCaloHit);
			}
		}

		float enteringPointX = 0.f;
		float enteringPointY = 0.f;

		for(unsigned int c=0 ; c<caloHitVecForEnteringPoint.size() ; c++)
		{
			CaloHit *pCaloHit = caloHitVecForEnteringPoint.at(c);

			enteringPointX += pCaloHit->GetIJK().at(0);
			enteringPointY += pCaloHit->GetIJK().at(1);
		}

		enteringPointX /= caloHitVecForEnteringPoint.size();
		enteringPointY /= caloHitVecForEnteringPoint.size();

		return ThreeVector(enteringPointX, enteringPointY, 0);
	}



}  // namespace 

