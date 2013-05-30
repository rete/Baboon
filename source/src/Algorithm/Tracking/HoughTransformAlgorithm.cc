/*
 *
 * HoughTransformAlgorithm.cc cpp file template generated by fclass
 * Creation date : Thu Mar 14 21:55:13 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#include "Algorithm/Tracking/HoughTransformAlgorithm.hh"


using namespace std;



namespace baboon {


	HoughTransformAlgorithm::HoughTransformAlgorithm()
		: AbstractAlgorithm("HoughTransformAlgorithm") {

		thetaMax = 100;
		rMax = 140;
		houghSpaceX = 0;
		houghSpaceY = 0;
		needData = true;
		needParams = true;

	}

	HoughTransformAlgorithm::~HoughTransformAlgorithm() {

		this->DeleteHoughSpace();
	}


	void HoughTransformAlgorithm::Init() {

		data.GetValue("thetaMax",&thetaMax);
		data.GetValue("rMax",&rMax);
		data.GetValue("clusterSizeLimit",&clusterSizeLimit);
		data.GetValue("minimumBinning",&minimumBinning);
		data.GetValue("deltaPosMax",&deltaPosMax);
		data.GetValue("trackSegmentMinimumSize",&trackSegmentMinimumSize);
		data.GetValue("maximumDistanceBetweenHitsInPlane",&maximumDistanceBetweenHitsInPlane);
		data.GetValue("maximumDistanceBetweenHitsForLayers",&maximumDistanceBetweenHitsForLayers);
		data.GetValue("maximumDistanceAllowedForHitsInTrack",&maximumDistanceAllowedForHitsInTrack);


		this->DeleteHoughSpace();
		this->AllocateHoughSpace();

	}

	void HoughTransformAlgorithm::AllocateHoughSpace() {

		houghSpaceX = new int *[thetaMax];
		for( int i=0 ;i<thetaMax ; i++ ) {
			houghSpaceX[i] = new int [rMax];
			for( int j=0 ;j<thetaMax ; j++ ) {
				houghSpaceX[i][j] = 0;
			}
		}

		houghSpaceY = new int *[thetaMax];
		for( int i=0 ;i<thetaMax ; i++ ) {
			houghSpaceY[i] = new int [rMax];
			for( int j=0 ;j<thetaMax ; j++ ) {
				houghSpaceY[i][j] = 0;
			}
		}
	}

	void HoughTransformAlgorithm::DeleteHoughSpace() {

		if( houghSpaceX != 0 ) {
			for(int i=0 ; i < thetaMax ; ++i) {
				delete [] houghSpaceX[i];
			}
			delete [] houghSpaceX;
		}

		if( houghSpaceY != 0 ) {
			for(int i=0 ; i < thetaMax ; ++i) {
				delete [] houghSpaceY[i];
			}
			delete [] houghSpaceY;
		}
	}


	Return HoughTransformAlgorithm::CheckConsistency() {

		return BABOON_SUCCESS();
	}


/////////////////////////////////////////////////////////////////
//			 _____________________________________________
//			< Beginning of the main part of the algorithm >
//			 ---------------------------------------------
//				 \
//				  \
//				   ("`-'  '-/") .___..--' ' "`-._
//					 ` *_ *  )    `-.   (      ) .`-.__. `)
//					 (_Y_.) ' ._   )   `._` ;  `` -. .-'
//				  _.. `--'_..-_/   /--' _ .' ,4
//			   ( i l ),-''  ( l i),'  ( ( ! .-'
////////////////////////////////////////////////////////////////


	void HoughTransformAlgorithm::Execute() {

		ClusterCollection *clusterCollection = ClusteringManager::GetInstance()->GetCluster2D();

		if( clusterCollection->empty() )
			return;

		HoughClusterCollection *houghClusterCollection = new HoughClusterCollection();

		for( unsigned int i=0 ; i<clusterCollection->size() ; i++ ) {

			Cluster *cluster = clusterCollection->at(i);
			if( cluster->GetClusterSize() > clusterSizeLimit ) continue;
			if( cluster->IsIsolatedFromClusters(clusterCollection) ) {

				HoughCluster *houghCluster = new HoughCluster();

				for( int t=0 ; t<thetaMax ; t++ ) {

					houghCluster->rhox.push_back( (int) abs( cluster->GetPosition().z() * cos(-M_PI/2.0 + t*M_PI/thetaMax)
									+ cluster->GetPosition().x() * sin(-M_PI/2.0 + t*M_PI/thetaMax) ) );
					houghCluster->rhoy.push_back( (int) abs( cluster->GetPosition().z() * cos(-M_PI/2.0 + t*M_PI/thetaMax)
									+ cluster->GetPosition().y() * sin(-M_PI/2.0 + t*M_PI/thetaMax) ) );
					houghSpaceX[ t ][ houghCluster->rhox.at( t ) ]++;
					
				}
				houghCluster->cluster = cluster;
				houghClusterCollection->push_back( houghCluster );
			}
		}

		for( int t=0 ; t<thetaMax ; t++ ) {
			for( unsigned int i=0 ; i<houghClusterCollection->size() ; i++ ) {

				HoughCluster *houghCluster = houghClusterCollection->at(i);
				if( houghSpaceX [ t ][ houghCluster->rhox.at(t) ] < minimumBinning ) continue;

				houghCluster->tagx = fGood;
				int inc = 0;

				for( unsigned int j=0 ; j<houghClusterCollection->size() ; j++ ) {

					HoughCluster *houghCluster2 = houghClusterCollection->at(j);
					if( houghCluster->cluster->GetPosition().z() == houghCluster2->cluster->GetPosition().z()
					 || houghCluster->rhox.at(t) - houghCluster2->rhox.at(t) != 0 )
						continue;

					if( abs( houghCluster2->cluster->GetPosition().x() - houghCluster->cluster->GetPosition().x() ) < deltaPosMax
					 && abs( houghCluster2->cluster->GetPosition().y() - houghCluster->cluster->GetPosition().y() ) < deltaPosMax
					 && abs( houghCluster2->cluster->GetPosition().z() - houghCluster->cluster->GetPosition().z() ) < deltaPosMax )
						inc++;
				}

				if( inc < 2 ) continue;

				for( int t2=0 ; t2<thetaMax ; t2++ )
					houghSpaceY[ t2 ][ houghCluster->rhoy.at(t2) ] ++;

			}
		}

		for( int t=0 ; t<thetaMax ; t++ ) {
			for( unsigned int i=0 ; i<houghClusterCollection->size() ; i++ ) {

				HoughCluster *houghCluster = houghClusterCollection->at(i);
				if( houghSpaceY [ t ][ houghCluster->rhoy.at(t) ] < minimumBinning ) continue;
				houghCluster->tagy = fGood;
			}
		}

		for( unsigned int i=0 ; i<houghClusterCollection->size() ; i++ ) {
			HoughCluster *houghCluster = houghClusterCollection->at(i);
			if( houghCluster->tagx == fGood && houghCluster->tagy == fGood ) {
				houghCluster->finalTag = fGood;
			}
		}

		TrackCollection *trackCollection = new TrackCollection();
		HoughClusterCollection houghClusterTemp;

		for( unsigned int i=0 ; i<houghClusterCollection->size() ; i++ ) {

			HoughCluster *houghCluster = houghClusterCollection->at( i );
			if( houghCluster->finalTag != fGood ) continue;

			if( find(houghClusterTemp.begin()
									,houghClusterTemp.end()
									,houghCluster )
									!= houghClusterTemp.end() )
								continue;

			houghClusterTemp.push_back( houghClusterCollection->at( i ) );

			HoughClusterCollection tracks;
			tracks.push_back( houghClusterCollection->at( i ) );

			for(unsigned int j=0 ; j<houghClusterCollection->size() ; j++) {

				HoughCluster *houghCluster2 = houghClusterCollection->at( j );

				if( houghCluster2->finalTag != fGood ) continue;

				if( find(houghClusterTemp.begin()
						,houghClusterTemp.end()
						,houghCluster2 )
						!= houghClusterTemp.end() )
					continue;

				for( unsigned int k=0 ; k<tracks.size() ; k++ ) {

					HoughCluster *track = tracks.at(k);

					if( abs(track->cluster->GetPosition().z() - houghCluster2->cluster->GetPosition().z()) <= maximumDistanceBetweenHitsForLayers
					 && abs(track->cluster->GetPosition().y() - houghCluster2->cluster->GetPosition().y()) <= maximumDistanceBetweenHitsInPlane
					 && abs(track->cluster->GetPosition().x() - houghCluster2->cluster->GetPosition().x()) <= maximumDistanceBetweenHitsInPlane ) {

						tracks.push_back( houghClusterCollection->at( j ) );
						houghClusterTemp.push_back( houghClusterCollection->at( j ) );
						break;
					}
				}
			}


			HitCollection tempHitCollection;

			// Copy all the clusters hits in a hit collection.
			for(unsigned int j=0 ; j<tracks.size() ; j++) {
				HitCollection *hitColTemp = tracks.at(j)->cluster->GetHitCollection();
				for(unsigned int k=0 ; k<hitColTemp->size() ; k++ ) {
					tempHitCollection.push_back( hitColTemp->at(k) );
				}
			}

			// Check for hits in track which are to far from the other hits in the track and remove them.
			for( unsigned int k=0 ; k<tempHitCollection.size() ; k++ ) {

				Hit *hit1 = tempHitCollection.at(k);
				IntVector ijk1 = hit1->GetIJK();
				double minimumDistance = 10000;

				for( unsigned int j=0 ; j<tempHitCollection.size() ; j++ ) {

					if( k == j ) continue;
					Hit *hit2 = tempHitCollection.at(j);
					IntVector ijk2 = hit2->GetIJK();

					double distance = sqrt( (ijk2.at(0)-ijk1.at(0)) * (ijk2.at(0)-ijk1.at(0))
								  + (ijk2.at(1)-ijk1.at(1)) * (ijk2.at(1)-ijk1.at(1))
								  + (ijk2.at(2)-ijk1.at(2)) * (ijk2.at(2)-ijk1.at(2)) );

					if( distance < minimumDistance )
						minimumDistance = distance;
				}
				if( minimumDistance > maximumDistanceAllowedForHitsInTrack ) {
					tempHitCollection.erase( tempHitCollection.begin() + k );
					k--;
				}
			}

			// Un-tag all the clusters of the track. No need anymore.
			for( unsigned int k=0 ; k<tracks.size() ; k++ )
				tracks.at(k)->cluster->SetClusterTag( fUndefined );

			/*
			 * Check if the number of remaining hits in the hit track
			 * collection is enough to build a track
			 */
			if( tempHitCollection.size() < trackSegmentMinimumSize ) {
				tempHitCollection.clear();
				continue;
			}

			// Build the track and tag the hits inside as "track hits"
			Track *track = new Track();
			for(unsigned int j=0 ; j<tempHitCollection.size() ; j++) {
				track->AddHit( tempHitCollection.at(j) );

			}


//			trackCollection->push_back( track );
			track->SortHits();

			if( abs(track->GetHitCollection()->at(0)->GetIJK().at(2) - track->GetHitCollection()->at( track->Size() -1 )->GetIJK().at(2) ) < 4 ) {

				delete track;
				continue;
			}

			for(unsigned int j=0 ; j<tempHitCollection.size() ; j++) {
				tempHitCollection.at(j)->SetHitTag( fTrack );
			}

			track->GetExtremities().first->SetHitTag( fTrackExtremity );
			track->GetExtremities().second->SetHitTag( fTrackExtremity );

//			vector<ThreeVector> positions = track->GetPositions();
//			vector<ThreeVector> weigths( track->Size() , ThreeVector(1,1,1) );

//			Linear3DFit fitter(positions,weigths);
//			fitter.Fit();
//			double chi2 = fitter.GetChi2();
//			cout << "chi2 : " << chi2 << endl;

			TrackManager::GetInstance()->AddTrack( track );
		}

		houghClusterTemp.clear();

		for( unsigned int i=0 ; i<houghClusterCollection->size() ; i++ )
			delete houghClusterCollection->at( i );

		houghClusterCollection->clear();
		delete houghClusterCollection;
	}



	void HoughTransformAlgorithm::End() {}

}


