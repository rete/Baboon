/*
 *
 * Cluster.cc cpp file template generated by fclass
 * Creation date : Sat Mar  2 02:06:41 2013
 * Copyright (c) 2012 CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#include "Objects/Cluster.hh"

using namespace std ;

namespace baboon {

	Cluster::Cluster()
		: HitCompositeObject(),
		  TypedObject("Cluster") {

		fClusterTag = UndefinedTag();
	}

	Cluster::~Cluster() {

	}


	void Cluster::ComputePosition( PositionComputation computation ) {

		if( caloHitCollection->empty() )
			return;

		ThreeVector pos(0,0,0);

		if( computation == fComputeCell ) {
			for(unsigned int i=0 ; i<caloHitCollection->size() ; i++) {

				CaloHit *caloHit = caloHitCollection->at(i);
				IntVector hitIJK = caloHit->GetIJK();

				pos.set( pos.x() + hitIJK.at(0)
						,pos.y() + hitIJK.at(1)
						,pos.z() + hitIJK.at(2) );
			}
		}
		else {
			for(unsigned int i=0 ; i<caloHitCollection->size() ; i++) {

				CaloHit *caloHit = caloHitCollection->at(i);
				ThreeVector hitPos = caloHit->GetPosition();

				pos.set( pos.x() + hitPos.x()
						,pos.y() + hitPos.y()
						,pos.z() + hitPos.z() );
			}
		}


		pos *= double(1.0/caloHitCollection->size());
		position = pos;
	}

	void Cluster::SetPosition( const ThreeVector &pos ) {

		position = pos;
	}

	bool Cluster::IsIsolatedFromClusters(const ClusterCollection* clusters) {

		bool isol=false;
		int neighbour=0;
		int neighbourbis=0;
		int neighbourbisbis=0;

		ComputePosition();
		for(int clustID=0 ; clustID<clusters->size() ; clustID++) {

			ThreeVector pos(clusters->at(clustID)->GetPosition() );

			if(abs(position.x()-pos.x())>1
			&& abs(position.x()-pos.x())<5
			&& abs(position.y()-pos.y())>1
			&& abs(position.y()-pos.y())<5
			&& position.z()==pos.z())
				{ neighbour++; }

			if(position.z()!=pos.z()
			&& abs(position.z()-pos.z())<3
			&& abs(position.x()-pos.x())<5
			&& abs(position.y()-pos.y())<5)
				{ neighbourbis++; }

			if(position.z()!=pos.z()
			&& abs(position.z()-pos.z())<3
			&& abs(position.x()-pos.x())<5
			&& abs(position.y()-pos.y())<5
			&& clusters->at(clustID)->Size()>8 )
				{ neighbourbisbis++; }
		}

		if( neighbour<2
		&& neighbourbis>0
		&& neighbourbis<6
		&& neighbourbisbis<2 ) {

			isol=true;
			}

		return isol;
	}

	void Cluster::SetClusterTag( const BaseTag &clustTag ) {

		fClusterTag.CopyTag( clustTag );
	}


	void Cluster::SetClusterTagRecursive( const BaseTag &clustTag ) {

		for(unsigned int i=0 ; i<caloHitCollection->size() ; i++) {
			caloHitCollection->at(i)->SetTag( clustTag );
		}
		this->SetClusterTag( clustTag );
	}


	Return Cluster::SetClusterType( const ClusterType &type ) {

		fType = type;
		return BABOON_SUCCESS();
	}



}
