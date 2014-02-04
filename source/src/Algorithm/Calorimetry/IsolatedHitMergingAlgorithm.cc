  /// \file IsolatedHitMergingAlgorithm.cc
/*
 *
 * IsolatedHitMergingAlgorithm.cc source template generated by fclass
 * Creation date : mer. d�c. 11 2013
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
 * @author : Remi Ete
 * @version
 * @copyright
 *
 *
 */


#include "Algorithm/Calorimetry/IsolatedHitMergingAlgorithm.hh"

namespace baboon {

	IsolatedHitMergingAlgorithm::IsolatedHitMergingAlgorithm()
		: AbstractAlgorithm("IsolatedHitMergingAlgorithm") {

		needData = false;
		calorimeter = 0;
	}

	IsolatedHitMergingAlgorithm::~IsolatedHitMergingAlgorithm() {

	}


	Return IsolatedHitMergingAlgorithm::Init() {

		return BABOON_SUCCESS();
	}


	Return IsolatedHitMergingAlgorithm::CheckConsistency() {

		BABOON_CHECK_POINTER( calorimeter );
		return BABOON_SUCCESS();
	}


	Return IsolatedHitMergingAlgorithm::Execute() {

		ClusterCollection *clusters = ClusteringManager::GetInstance()->GetCluster3D();

		if( clusters->empty() )
			return BABOON_SUCCESS();

		CaloHitCollection *caloHitCollection = calorimeter->GetCaloHitCollection();

		if( caloHitCollection->empty() )
			return BABOON_SUCCESS();

		std::map<CaloHit*,Cluster*> caloHitToClusterAssociation;

		for( unsigned int c1=0 ; c1<caloHitCollection->size() ; c1++ ) {

			CaloHit *caloHit = caloHitCollection->at(c1);

			if( caloHit->GetTag() != IsolatedTag() )
				continue;

			double minimumDistance;
			Cluster *closestCluster = 0;

			for( unsigned int cl=0 ; cl<clusters->size() ; cl++ ) {

				if( closestCluster == 0 ) {

					closestCluster = clusters->at(cl);
					minimumDistance = DistanceToCluster( caloHit->GetPosition() , clusters->at(cl) );
//					minimumDistance *= 1 + std::log (((double)biggestClusterSize/(double)clusters->at(cl)->Size()) );
					continue;
				}

				double distance = DistanceToCluster( caloHit->GetPosition() , clusters->at(cl) );
//				distance *= 1 + std::log( ((double)biggestClusterSize/(double)clusters->at(cl)->Size()) );

				if( distance < minimumDistance ) {
					minimumDistance = distance;
					closestCluster = clusters->at(cl);
				}
			}

			caloHitToClusterAssociation[ caloHit ] = closestCluster;
		}

		for( std::map<CaloHit*,Cluster*>::iterator it=caloHitToClusterAssociation.begin() ;
			it!=caloHitToClusterAssociation.end() ; it++ ) {

			it->second->AddCaloHit( it->first );
		}
		caloHitToClusterAssociation.clear();


		// for monitoring
		for( unsigned int cl=0 ; cl<clusters->size() ; cl++ ) {

			CaloHitCollection *clusterHits = clusters->at(cl)->GetCaloHitCollection();

			for( unsigned int c1=0 ; c1<clusterHits->size() ; c1++ ) {
				clusterHits->at(c1)->SetColor( cl + 1 );
			}
		}

		return BABOON_SUCCESS();
	}


	Return IsolatedHitMergingAlgorithm::End() {

		calorimeter = 0;
		return BABOON_SUCCESS();
	}

// ------------------------------------------------------------------------------------------------------------------------

	double IsolatedHitMergingAlgorithm::DistanceToCluster( const ThreeVector &pos , Cluster *cluster ) {

		if( cluster == 0 )
			return 0.0;

		CaloHitCollection *clusterHits = cluster->GetCaloHitCollection();

		double distanceMinimum = 10000000.0;

		for( unsigned int i=0 ; i<clusterHits->size() ; i++ ) {

			double distance = (pos - clusterHits->at(i)->GetPosition() ).mag();// / (double)cluster->Size();

			if( distance < distanceMinimum )
				distanceMinimum = distance;

		}

		return distanceMinimum;
	}


}  // namespace 

