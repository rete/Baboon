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
	}

	ClusteringAlgorithm::~ClusteringAlgorithm() {

	}


	void ClusteringAlgorithm::Init() {

	 }


	void ClusteringAlgorithm::Execute() {

		HitManager *hitManager = HitManager::GetInstance();
		HitCollection *hitCollection = hitManager->GetHitCollection();
		HitCollection treatedHits;

		for( unsigned int hitID=0 ; hitID<hitCollection->size() ; hitID++ ) {

			if( find( treatedHits.begin() , treatedHits.end() , hitCollection->at(hitID) ) != treatedHits.end() )
				continue;

			if ( !CheckTag( hitCollection->at(hitID)->GetHitTag() ) ) continue;

			treatedHits.push_back( hitCollection->at(hitID) );
			HitCollection *hitCol = new HitCollection();
			hitCol->push_back( hitCollection->at(hitID) );
			IntVector ijk = hitCollection->at(hitID)->GetIJK();

			for( unsigned int hitID2=0 ; hitID2<hitCollection->size() ; hitID2++ ) {

				Hit *hit2 = hitCollection->at(hitID2);

				if ( !CheckTag( hit2->GetHitTag() ) )
					continue;

				if( find( treatedHits.begin() , treatedHits.end() , hit2 ) != treatedHits.end() )
					continue;

				IntVector ijk2 = hitCollection->at(hitID2)->GetIJK();

				if( fClusteringMode == fClustering2D )
					if( abs( ijk.at(2) == ijk2.at(2) ) )
						continue;


				for( unsigned int hitID3=0 ; hitID3<hitCol->size() ; hitID3++ ) {

					IntVector ijk3 = hitCol->at(hitID3)->GetIJK();

					if( abs( ijk3.at(0)-ijk2.at(0) ) < 2
					 && abs( ijk3.at(1)-ijk2.at(1) ) < 2
					 && abs( ijk3.at(2)-ijk2.at(2) ) < 2 ) {

						treatedHits.push_back( hit2 );
						hitCol->push_back( hit2 );
						break;
					}
				}
			}

			Cluster *cluster = new Cluster();
			cluster->SetHitCollection( hitCol );
			clusterCollection->push_back( cluster );
//			getchar();
		}
		treatedHits.clear();

	}


	void ClusteringAlgorithm::End() {

		hitTagToCluster.clear();
		hitTagToAvoid.clear();
	}


	Return ClusteringAlgorithm::CheckConsistency() {

		Return ret;
		if( clusterCollection == 0 || clusterCollection == NULL )
			return S_ERROR("While checking consistency : Cluster collection not set or set to 0 in ClusteringAlgorithm. ");
		if( !clusterCollection->empty() ) {
			for( unsigned int i=0 ; i<clusterCollection->size() ; i++ )
				delete clusterCollection->at(i);
			clusterCollection->clear();
			ret = S_OK("Warning : ClusterCollection has been cleared while checking consistency");
		}
		for( unsigned int i=0 ; i<hitTagToCluster.size() ; i++ ) {
			if( std::find( hitTagToAvoid.begin() , hitTagToAvoid.end() , hitTagToCluster.at(i) ) != hitTagToAvoid.end() ) {
				return S_ERROR("Tag to be avoided also set to be clustered ... Check your inputs!");
			}
		}
		return ret;
	}


	Return ClusteringAlgorithm::AddHitTagToCluster( const Tag &fTag ) {

		if( std::find( hitTagToCluster.begin() , hitTagToCluster.end() , fTag ) == hitTagToCluster.end() ) {
			hitTagToCluster.push_back( fTag );
			return S_OK();
		}
		else return S_OK("Warning : hit tag was already added");
	}


	Return ClusteringAlgorithm::AddHitTagToAvoid( const Tag &fTag ) {

		if( std::find( hitTagToAvoid.begin() , hitTagToAvoid.end() , fTag ) == hitTagToAvoid.end() ) {
			hitTagToAvoid.push_back( fTag );
			return S_OK();
		}
		else return S_OK("Warning : hit tag was already added");
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


}  // namespace 

