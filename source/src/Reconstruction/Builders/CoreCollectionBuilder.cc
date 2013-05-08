  /// \file CoreCollectionBuilder.cc
/*
 *
 * CoreCollectionBuilder.cc source template generated by fclass
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


#include "Reconstruction/Builders/CoreCollectionBuilder.hh"


using namespace std;


namespace baboon {


	CoreCollectionBuilder *CoreCollectionBuilder::instance = 0;


	CoreCollectionBuilder::CoreCollectionBuilder()
		: ObjectBuilder<CoreCollection>() {

	}


	CoreCollectionBuilder::~CoreCollectionBuilder() {

	}

	CoreCollectionBuilder *CoreCollectionBuilder::GetInstance() {

		if( instance == 0 ) instance = new CoreCollectionBuilder();
		return instance;
	}

	void CoreCollectionBuilder::Kill() {

		if( instance != 0 ) {
			delete instance;
			instance = 0;
		}
	}

	Return CoreCollectionBuilder::ClearObject() {

		if( object != 0 ) {
			if( object->empty() ) {
				delete object;
				object = 0;
				return S_OK();
			}

			for( unsigned int i=0 ; i<object->size() ; i++ ) {
				if(object->at(i) != 0) delete object->at(i);
				object->at(i) = 0;
			}
			object->clear();
			delete object;
			object = 0;
		}

		return S_OK();
	}

	void CoreCollectionBuilder::BuildObject() {

		if( object != 0) ClearObject();

		object = new CoreCollection();

		ClusterCollection *clusterCollection = new ClusterCollection();

		ClusteringAlgorithm *clustAlgo = new ClusteringAlgorithm();
		clustAlgo->SetClusterCollection( clusterCollection );
		clustAlgo->SetHitTagToCluster( fCore );
		clustAlgo->Process();
		delete clustAlgo;

		for( unsigned int i=0 ; i<clusterCollection->size() ; i++ ) {

			Core *core = new Core();
			HitCollection *hitCollection = clusterCollection->at(i)->GetHitCollection();

			for( unsigned int j=0 ; j<hitCollection->size() ; j++ )
				core->AddHit( hitCollection->at(j) );

			delete clusterCollection->at(i);
			object->push_back( core );
		}
		clusterCollection->clear();
		delete clusterCollection;


		cout << "number of cores : " << object->size() << endl;

	}



	// template initialization
	template ObjectBuilder<CoreCollection>::ObjectBuilder();
	template Return ObjectBuilder<CoreCollection>::SetObject( CoreCollection *obj );


}  // namespace 
