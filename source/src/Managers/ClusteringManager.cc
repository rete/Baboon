/*
 *
 * ClusteringManager.cc cpp file template generated by fclass
 * Creation date : Fri Mar 15 18:06:44 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#include "Managers/ClusteringManager.hh"



using namespace std ;


namespace baboon {

	ClusteringManager *ClusteringManager::instance = 0;

	ClusteringManager::ClusteringManager() {

		clusters3D = new ClusterCollection();
		clusters2D = new ClusterCollection();
	}

	ClusteringManager::~ClusteringManager() {

		delete clusters3D;
		delete clusters2D;
	}


	ClusteringManager *ClusteringManager::GetInstance() {

		if(instance == 0)
			instance = new ClusteringManager;
		return instance;
	}

	void ClusteringManager::Kill() {
		if(instance != NULL) {
			delete instance;
			instance = NULL;
		}
	}


	Return ClusteringManager::AddCluster( Cluster *cluster ) {

		if( cluster == 0 )
			return BABOON_INVALID_PARAMETER("While adding a cluster : assertion cluster != 0 failed");

		if( cluster->GetClusterType() == fCluster2D ) {
			ClusterCollection::iterator clusterIt = std::find( clusters2D->begin() ,clusters2D->end() , cluster );
			if( clusterIt != clusters2D->end() )
				return BABOON_ALREADY_PRESENT("Cluster already exists in cluster collection");
			else {
				clusters2D->push_back( cluster );
				return BABOON_SUCCESS();
			}
		}
		else if( cluster->GetClusterType() == fCluster3D ) {
			ClusterCollection::iterator clusterIt = std::find( clusters3D->begin() ,clusters3D->end() , cluster );
			if( clusterIt != clusters3D->end() )
				return BABOON_ALREADY_PRESENT("Cluster already exists in cluster collection");
			else {
				clusters3D->push_back( cluster );
				return BABOON_SUCCESS();
			}
		}
		return BABOON_ERROR("Cluster type undefined...");
	}


	Return ClusteringManager::RemoveCluster( Cluster *cluster ) {

		if( cluster == 0 )
			return BABOON_INVALID_PARAMETER("While removing a cluster : assertion cluster != 0 failed");

		if( cluster->GetClusterType() == fCluster2D ) {
			ClusterCollection::iterator clusterIt = std::find( clusters2D->begin() ,clusters2D->end() , cluster );
			if( clusterIt != clusters2D->end() ) {
				delete cluster;
				clusters2D->erase(clusterIt);
			}
			else return BABOON_NOT_FOUND("While removing cluster : cluster was not registered in the cluster collection 2D");
		}
		else if( cluster->GetClusterType() == fCluster3D ) {
			ClusterCollection::iterator clusterIt = std::find( clusters3D->begin() ,clusters3D->end() , cluster );
			if( clusterIt != clusters3D->end() ) {
				delete cluster;
				clusters3D->erase(clusterIt);
			}
			else return BABOON_NOT_FOUND("While removing cluster : cluster was not registered in the cluster collection 3D");
		}
		return BABOON_ERROR("Cluster type undefined...");
	}


	Return ClusteringManager::ClearAllContent() {

		if( clusters2D == 0 )
			return BABOON_INVALID_PARAMETER("While clearing all content in clustering manager : assertion clusters2D != 0 failed");
		if( clusters3D == 0 )
			return BABOON_INVALID_PARAMETER("While clearing all content in clustering manager : assertion clusters3D != 0 failed");

		for( unsigned int i=0 ; i<clusters2D->size() ; i++ ) {
			if( clusters2D->at(i) != 0 )
				delete clusters2D->at(i);
		}
		for( unsigned int i=0 ; i<clusters3D->size() ; i++ ) {
			if( clusters3D->at(i) != 0 )
				delete clusters3D->at(i);
		}
		clusters2D->clear();
		clusters3D->clear();

		return BABOON_SUCCESS("Content cleared in clustering manager");
	}


	bool ClusteringManager::ClusterContainsHit( Cluster *cluster , CaloHit *caloHit ) {

		return cluster->Contains( caloHit );
	}


	Return ClusteringManager::MergeAndDeleteClusters( Cluster *clusterToEnlarge , Cluster *clusterToDelete ) {

		if( clusterToEnlarge == 0 || clusterToDelete == 0 )
			return BABOON_INVALID_PARAMETER("Assertion cluster != 0 failed");

		if( clusterToEnlarge == clusterToDelete )
			return BABOON_INVALID_PARAMETER("Cluster are the same. Can't merge the same clusters!");

		CaloHitCollection *caloHitCollection = clusterToDelete->GetCaloHitCollection();

		for( unsigned int i=0 ; i<caloHitCollection->size() ; i++ ) {

			clusterToEnlarge->AddCaloHit( caloHitCollection->at(i) );
			clusterToDelete->RemoveCaloHit( caloHitCollection->at(i) );
		}

		BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , this->RemoveCluster( clusterToDelete ) );

		delete clusterToDelete;

		return BABOON_SUCCESS();
	}


	Cluster *ClusteringManager::GetClusterAt( ClusterType clusterType , unsigned int I , unsigned int J , unsigned int K ) {

		CaloHit *CaloHit = 0;
//
//		if( clusterType == fCluster2D ) {
//			for( unsigned int i=0 ; i<clusters2D->size() ; i++ )
//				if( clusters2D->at(i)->Contains( hitAtIJK ) )
//					return clusters2D->at(i);
//		}
//		else {
//			for( unsigned int i=0 ; i<clusters3D->size() ; i++ )
//				if( clusters3D->at(i)->Contains( hitAtIJK ) )
//					return clusters3D->at(i);
//		}

		return 0;
	}


	Cluster *ClusteringManager::GetClusterContainingCaloHit( ClusterType clusterType , CaloHit *caloHit ) {

		if( caloHit == 0 )
			return 0;

		if( clusterType == fCluster2D ) {
			for( unsigned int i=0 ; i<clusters2D->size() ; i++ )
				if( clusters2D->at(i)->Contains( caloHit ) )
					return clusters2D->at(i);
		}
		else {
			for( unsigned int i=0 ; i<clusters3D->size() ; i++ )
				if( clusters3D->at(i)->Contains( caloHit ) )
					return clusters3D->at(i);
		}

		return 0;
	}


}





