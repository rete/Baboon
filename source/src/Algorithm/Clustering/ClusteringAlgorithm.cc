  /// \file ClusteringAlgorithm.cc
/*
 *
 * ClusteringAlgorithm.cc source template generated by fclass
 * Creation date : mar. mai 7 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Algorithm/Clustering/ClusteringAlgorithm.hh"

using namespace std;


namespace baboon {



	ClusteringAlgorithm::ClusteringAlgorithm()
		: AbstractAlgorithm("ClusteringAlgorithm") {

		needData = false;
		fClusteringMode = fClustering2D;
		fTaggingMode = fClusterTagMode;
		clusterCollection = 0;
		clusterSizeLowerLimit = 1;
		neighborDistance = 1;

	}

	ClusteringAlgorithm::~ClusteringAlgorithm() {

	}


	Return ClusteringAlgorithm::Init() {


		return BABOON_SUCCESS();
	 }


	int IJKToKey( int I , int J , int K ) {
		return 100*100*K+100*J+I;
	}
	std::vector<int> decoder(int ijk){
		std::vector<int> dec;
		dec.push_back(ijk%100);
		dec.push_back((ijk/100)%100);
		dec.push_back(ijk/100/100);
		return dec;
	}

	Return ClusteringAlgorithm::Execute() {

		HitManager *hitManager = HitManager::GetInstance();
		HitCollection *hitCollection = hitManager->GetHitCollection();

		for( unsigned int i=0 ; i<hitCollection->size() ; i++ ) {

			Hit *hit = hitCollection->at(i);
			IntVector ijk = hit->GetIJK();

			if( find( treatedHits.begin() , treatedHits.end() , hit ) != treatedHits.end() )
				continue;

//			cout << "Hit at : " << ijk.at(0) << " " << ijk.at(1) << " " << ijk.at(2) << endl;

			treatedHits.push_back( hit );

			if ( !CheckTag( hitCollection->at( i )->GetHitTag() ) )
				continue;

			Cluster *cluster = new Cluster();
			cluster->AddHit( hit );

			this->FindCluster( hit , cluster );

			if( fClusteringMode == fClustering2D ) cluster->SetType( fCluster2D );
			else if( fClusteringMode == fClustering3D ) cluster->SetType( fCluster3D );

			if( cluster->GetClusterSize() < clusterSizeLowerLimit ) {
				delete cluster;
				continue;
			}

			clusterCollection->push_back( cluster );

		}

		/*
		map<int,Hit*> hitMap;

		for( unsigned int i=0 ; i<hitCollection->size(); i++ ) {

			IntVector ijk = hitCollection->at(i)->GetIJK();

			int key = IJKToKey( ijk.at(0), ijk.at(1) , ijk.at(2)  );
			hitMap[  key ] = hitCollection->at(i);
		}
		//for(std::map<int,Hit*>::iterator it=hitMap.begin(); it!=hitMap.end(); ++it)
			//std::cout << it->first <<"=?"<< it->second->GetIJK()[2]<<it->second->GetIJK()[1]<<it->second->GetIJK()[0]<< std::endl;


		HitCollection treatedHits;

		map<int,Hit*>::iterator it1;
//		for( unsigned int hitID=0 ; hitID<hitCollection->size() ; hitID++ ) {

		for( it1=hitMap.begin() ; it1!=hitMap.end() ; it1++ ) {

			if( find( treatedHits.begin() , treatedHits.end() , it1->second ) != treatedHits.end() )
				continue;

//			if ( !CheckTag( hitCollection->at(hitID)->GetHitTag() ) ) continue;

			treatedHits.push_back( it1->second );
			HitCollection *hitCol = new HitCollection();
			hitCol->push_back( it1->second );
			IntVector ijk = it1->second->GetIJK();
			//std::vector<int> ijk=decoder(it1>first);
			map<int,Hit*>::iterator it2;
			for( it2=hitMap.begin() ; it2!=hitMap.end() ; it2++ ) {

//				Hit *hit2 = hitCollection->at(hitID2);

//				if ( !CheckTag( hit2->GetHitTag() ) )
//					continue;

				if( find( treatedHits.begin() , treatedHits.end() , it2->second ) != treatedHits.end() )
					continue;

				IntVector ijk2 = it2->second->GetIJK();
				//std::vector<int> ijk2=decoder(it2->first);
//				if( fClusteringMode == fClustering2D )
					if(  ijk.at(2) != ijk2.at(2) )
						continue;


				for( unsigned int hitID3=0 ; hitID3<hitCol->size() ; hitID3++ ) {

					IntVector ijk3 = hitCol->at(hitID3)->GetIJK();
					if( ijk2.at(0)==49&&ijk2.at(1)==46&&ijk2.at(2)==2 ){
						std::cout << ijk3.at(0)<<", "<<ijk3.at(1)<<", "<<ijk3.at(2)<<std::endl;
						std::cout << ijk3.at(0)-ijk2.at(0)<<", "<<ijk3.at(1)-ijk2.at(1)<<", "<<ijk3.at(2)-ijk2.at(2)<<std::endl;
					}
					if( abs( ijk3.at(0)-ijk2.at(0) ) <= neighborDistance
					 && abs( ijk3.at(1)-ijk2.at(1) ) <= neighborDistance
//					 && abs( ijk3.at(2)-ijk2.at(2) ) <= neighborDistance
					 ) {

						treatedHits.push_back( it2->second );
						hitCol->push_back( it2->second );
						break;
					}
				}
			}
//			if( hitCol->size() < clusterSizeLowerLimit ) {
//				hitCol->clear();
//				delete hitCol;
//				continue;
//			}
			Cluster *cluster = new Cluster();
			cluster->SetHitCollection( hitCol );
			if( fClusteringMode == fClustering2D ) cluster->SetType( fCluster2D );
			else if( fClusteringMode == fClustering3D ) cluster->SetType( fCluster3D );
			clusterCollection->push_back( cluster );
		}
		treatedHits.clear();
		*/

		return BABOON_SUCCESS();
	}


	Return ClusteringAlgorithm::End() {

		treatedHits.clear();
		hitTagToCluster.clear();
		hitTagToAvoid.clear();
		neighborDistance = 1;
		clusterSizeLowerLimit = 1;
		return BABOON_SUCCESS();
	}

	Return ClusteringAlgorithm::FindCluster( Hit *hit , Cluster *cluster ) {

		HitManager *hitManager = HitManager::GetInstance();
		IntVector ijk = hit->GetIJK();

		int distance = neighborDistance;

		for( int i=-distance ; i<=distance ; i++ ) {
			for( int j=-distance ; j<=distance ; j++ ) {
				for( int k=-distance ; k<=distance ; k++ ) {

					if( fClusteringMode == fClustering2D )
						if( k != 0 )
							continue;

					if( !hitManager->PadExists( ijk.at(0)+i , ijk.at(1)+j , ijk.at(2)+k ) )
						continue;

					if( !hitManager->PadIsTouched( ijk.at(0)+i , ijk.at(1)+j , ijk.at(2)+k ) )
						continue;

					Hit *otherHit = hitManager->GetHitAt( ijk.at(0)+i , ijk.at(1)+j , ijk.at(2)+k );

					if( find( treatedHits.begin() , treatedHits.end() , otherHit ) != treatedHits.end() )
						continue;

					if( cluster->Contains( otherHit ) )
						continue;

					treatedHits.push_back( otherHit );

					if ( !CheckTag( otherHit->GetHitTag() ) )
						continue;

					cluster->AddHit( otherHit );
					this->FindCluster( otherHit , cluster );
				}
			}
		}
		return BABOON_SUCCESS();
	}

	Return ClusteringAlgorithm::CheckConsistency() {

		Return ret = BABOON_SUCCESS();
		if( clusterCollection == 0 || clusterCollection == NULL )
			return BABOON_ERROR("While checking consistency : Cluster collection not set or set to 0 in ClusteringAlgorithm. ");
		if( !clusterCollection->empty() ) {
			for( unsigned int i=0 ; i<clusterCollection->size() ; i++ )
				delete clusterCollection->at(i);
			clusterCollection->clear();
			ret = BABOON_SUCCESS("Warning : ClusterCollection has been cleared while checking consistency");
		}
		for( unsigned int i=0 ; i<hitTagToCluster.size() ; i++ ) {
			if( std::find( hitTagToAvoid.begin() , hitTagToAvoid.end() , hitTagToCluster.at(i) ) != hitTagToAvoid.end() ) {
				return BABOON_ERROR("Tag to be avoided also set to be clustered ... Check your inputs!");
			}
		}
		return ret;
	}


	Return ClusteringAlgorithm::AddHitTagToCluster( const Tag &fTag ) {

		if( std::find( hitTagToCluster.begin() , hitTagToCluster.end() , fTag ) == hitTagToCluster.end() ) {
			hitTagToCluster.push_back( fTag );
			return BABOON_SUCCESS();
		}
		return BABOON_ALREADY_PRESENT("Warning : hit tag was already added");
	}


	Return ClusteringAlgorithm::AddHitTagToAvoid( const Tag &fTag ) {

		if( std::find( hitTagToAvoid.begin() , hitTagToAvoid.end() , fTag ) == hitTagToAvoid.end() ) {
			hitTagToAvoid.push_back( fTag );
			return BABOON_SUCCESS();
		}
		else return BABOON_ALREADY_PRESENT("Warning : hit tag was already added");
	}


	bool ClusteringAlgorithm::AvoidTag( const Tag &fTag ) {

		if( std::find( hitTagToAvoid.begin() , hitTagToAvoid.end() , fTag ) != hitTagToAvoid.end() ) {
			return true;
		}
		return false;
	}


	bool ClusteringAlgorithm::KeepTag( const Tag &fTag ) {

		if( std::find( hitTagToCluster.begin() , hitTagToCluster.end() , fTag ) != hitTagToCluster.end() ) {
			return true;
		}
		return false;
	}

	bool ClusteringAlgorithm::CheckTag( const Tag &fTag ) {

		if( fTaggingMode == fAvoidTagMode ) {
			if( hitTagToAvoid.empty() ) return true;
			else return AvoidTag( fTag );
		}
		else if( fTaggingMode == fClusterTagMode ) {
			if( hitTagToCluster.empty() ) return true;
			else return KeepTag( fTag );
		}

	}

	Return ClusteringAlgorithm::SetClusterSizeLowerLimit( unsigned int limit ) {

		clusterSizeLowerLimit = limit;
		return BABOON_SUCCESS();
	}

	Return ClusteringAlgorithm::SetNeighborDistance( unsigned int distance ) {

		neighborDistance = distance;
		return BABOON_SUCCESS();
	}

	Return ClusteringAlgorithm::SetClusteringMode( ClusteringMode mode ) {

		fClusteringMode = mode;
		return BABOON_SUCCESS();
	}

	Return ClusteringAlgorithm::SetTaggingMode( TaggingMode mode ) {

		fTaggingMode = mode;
		return BABOON_SUCCESS();
	}

}  // namespace 

