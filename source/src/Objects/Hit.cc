/*
 *
 * Hit.cc cpp file template generated by fclass
 * Creation date : Wed Mar  6 17:01:13 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#include "Objects/Hit.hh"
#include "Objects/Cluster.hh"

using namespace std ;

namespace sdhcal {

	Hit::Hit() : IMPL::CalorimeterHitImpl() {

		hitTag = fUndefined;
		SdhcalConfig::GetInstance()->GetData("general").GetValue("codingPattern",&codingPattern);
		hitTag = fUndefined;
		fThreshold = fThresholdUndefined;
	    cluster3D = new Cluster();
	    cluster3D->AddHit(this);
	    cluster2D = new Cluster();
	    cluster2D->AddHit(this);
	}

	Hit::Hit(EVENT::CalorimeterHit* hit) {

		SdhcalConfig::GetInstance()->GetData("general").GetValue("codingPattern",&codingPattern);
		hitTag = fUndefined;
		if(hit == NULL) cout << "hit pointer null" << endl;
	    _cellID0 = hit->getCellID0();
	    _cellID1 = hit->getCellID1();
	    _energy = hit->getEnergy();
	    _energyError = hit->getEnergyError();
	    _time = hit->getTime();
	    position = ThreeVector( hit->getPosition()[0] , hit->getPosition()[1] , hit->getPosition()[2] );
	    _type = hit->getType();
//	    _rawHit = hit->getRawHit()->clone();
	    hitTag = fUndefined;
	    if(_energy == 2.0) fThreshold = fThreshold1;
	    if(_energy == 1.0) fThreshold = fThreshold2;
	    if(_energy == 3.0) fThreshold = fThreshold3;
	    cluster3D = new Cluster();
	    cluster3D->AddHit(this);
	    cluster2D = new Cluster();
	    cluster2D->AddHit(this);

	}

	Hit::~Hit() {
		delete cluster2D;
		delete cluster3D;
	}


	EVENT::IntVec Hit::GetIJK() {

		EVENT::IntVec vec;
		UTIL::CellIDDecoder<Hit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");
		vec.push_back( idDecoder(this)["I"] );
		vec.push_back( idDecoder(this)["J"] );
		vec.push_back( idDecoder(this)["K-1"] );
		return vec;
	}

	void Hit::SetIJK( int I , int J , int K ) {

		IMPL::LCCollectionVec *lcCol = new IMPL::LCCollectionVec( "HCALBarrel" );
		lcCol->addElement( this );
		UTIL::CellIDEncoder<IMPL::CalorimeterHitImpl> idEncoder(codingPattern,lcCol);
		idEncoder["I"] = I;
		idEncoder["J"] = J;
		idEncoder["K-1"] = K;
		idEncoder.setCellID( this );
		lcCol->removeElementAt(0);
		delete lcCol;

	}

	IMPL::CalorimeterHitImpl *Hit::ToCalorimeterHitImpl() {


		IMPL::CalorimeterHitImpl *calHit = new IMPL::CalorimeterHitImpl();
		calHit->setCellID0( this->getCellID0() );
		calHit->setCellID1( this->getCellID1() );
		calHit->setEnergy( this->getEnergy() );
		calHit->setEnergyError( this->getEnergyError() );
		calHit->setTime( this->getTime() );
		float pos[3];
		pos[0] = this->GetPosition().x();
		pos[1] = this->GetPosition().y();
		pos[2] = this->GetPosition().z();
		calHit->setPosition( pos );
		calHit->setType( this->getType() );
//		calHit->setRawHit( this->getRawHit()->clone() );

		return calHit;
	}



	bool Hit::IsIsolatedFromHits( const HitCollection* hitCol ) {

		unsigned int isolCond = 3;

		for(unsigned int i=0 ; i<hitCol->size() ; i++) {

			Hit *hit1 = hitCol->at(i);

			unsigned int inc = 0;
			while (inc != isolCond ) {



				inc++;
			}

		}


	}


	void Hit::SetCluster3D( Cluster *cl ) {
		cluster3D = cl;
	}

	void Hit::SetCluster2D( Cluster *cl ) {
		cluster2D = cl;
	}

	void Hit::MergeClusters3D( Hit *hit ) {

		HitCollection *hitCol = hit->GetCluster3D()->GetHitCollection();

		for( unsigned int i=0 ; i<hitCol->size() ; i++ ) {
			cluster3D->AddHit(hitCol->at(i));
		}
		if( cluster3D != hitCol->at(0)->GetCluster3D() ) {
			delete hitCol->at(0)->GetCluster3D();
			for( unsigned int i=0 ; i<hitCol->size() ; i++ ) hitCol->at(i)->SetCluster3D(cluster3D);
		}

	}

	void Hit::MergeClusters2D( Hit *hit ) {

		HitCollection *hitCol = hit->GetCluster2D()->GetHitCollection();

		for( unsigned int i=0 ; i<hitCol->size() ; i++ ) {
			cluster2D->AddHit(hitCol->at(i));
		}
		if( cluster2D != hitCol->at(0)->GetCluster2D() ) {
			delete hitCol->at(0)->GetCluster2D();
			for( unsigned int i=0 ; i<hitCol->size() ; i++ ) hitCol->at(i)->SetCluster2D(cluster2D);
		}

	}


}

