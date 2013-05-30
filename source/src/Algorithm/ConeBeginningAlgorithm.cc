  /// \file ConeBeginningAlgorithm.cc
/*
 *
 * ConeBeginningAlgorithm.cc source template generated by fclass
 * Creation date : ven. avr. 19 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Algorithm/ConeBeginningAlgorithm.hh"

using namespace std;

namespace baboon {

	ConeBeginningAlgorithm::ConeBeginningAlgorithm()
		: AbstractAlgorithm("ConeBeginningAlgorithm") {
		needData = true;
	}

	ConeBeginningAlgorithm::~ConeBeginningAlgorithm() {

	}



	void ConeBeginningAlgorithm::Init() {

		hitCollection = HitManager::GetInstance()->GetHitCollection();
		data.GetValue("coneMinimumOpeningAngle",&coneMinimumOpeningAngle);
		data.GetValue("coneMaximumOpeningAngle",&coneMaximumOpeningAngle);
		data.GetValue("angleStep",&angleStep);
		data.GetValue("coneLength",&coneLength);
		data.GetValue("coneDirectedVectorHalfRange",&coneDirectedVectorHalfRange);
		data.GetValue("coneDirectedVectorStep",&coneDirectedVectorStep);
		data.GetValue("conePadsLowerLimitCut",&conePadsLowerLimitCut);
		data.GetValue("coreSizeLowerLimit",&coreSizeLowerLimit);
		data.GetValue("coneBackwardDistance",&coneBackwardDistance);
		data.GetValue("coneStepDistance",&coneStepDistance);
		data.GetValue("clusterSizeAtPeakPositionUpperLimit",&clusterSizeAtPeakPositionUpperLimit);

		if( !nodeList.empty() ) {
			for( unsigned int i=0 ; i<nodeList.size() ; i++ ) {
				SDHCALPrototype::GetGeoManager()->GetTopVolume()->RemoveNode( nodeList.at(i) );
			}
			nodeList.clear();
		}

	}


	Return ConeBeginningAlgorithm::CheckConsistency() {

		if( hitCollection == 0 )
			return BABOON_ERROR("ConeBeginning Algorithm bad init. Please check your inputs!");
		return BABOON_SUCCESS();
	}


	void ConeBeginningAlgorithm::Execute() {

		HitManager* hitManager = HitManager::GetInstance();
		CoreManager *coreManager = CoreManager::GetInstance();
		ClusteringManager *clusteringManager = ClusteringManager::GetInstance();

		CoreCollection *coreCollection = coreManager->GetCoreCollection();
		if( coreCollection->empty() )
			return;

		ostringstream ss;

		for( unsigned int c=0 ; c<coreCollection->size() ; c++ ) {

			Core *core = coreCollection->at(c);

			if( core->Size() < coreSizeLowerLimit )
				continue;

			ThreeVector coreCog = core->CenterOfGravity();

			int firstHitLayer = core->GetFirstHitLayer();

			HitCollection *hitCol = core->GetHitCollection();
			DoubleVector coreRMS = hitCol->GetRMS();


			int maximumClusterSizeInCore = 0;
//			ClusterCollection clusterTempList;
			int maximumClusterSizeLayer = 0;
			HitCollection hitTempCollection;

			double maxRadius = 0;

			for( unsigned int i=0 ; i<hitCol->size() ; i++ ) {

				Hit *hit1 = hitCol->at(i);
				if( std::find( hitTempCollection.begin() , hitTempCollection.end() , hit1 ) != hitTempCollection.end() )
					continue;

				HitCollection corePlanCluster;
				hitTempCollection.push_back( hit1 );
				corePlanCluster.push_back( hit1 );
				int coreClusterSize = 0;

				for(  unsigned int j=0 ; j<hitCol->size() ; j++ ) {

					Hit *hit2 = hitCol->at(j);
					if( hit1 == hit2 ) {
						coreClusterSize++;
						continue;
					}
					if( std::find( hitTempCollection.begin() , hitTempCollection.end() , hit2 ) != hitTempCollection.end() )
						continue;

					if( hit1->GetIJK().at(2) == hit2->GetIJK().at(2) ) {
						corePlanCluster.push_back( hit2 );
						coreClusterSize++;
						hitTempCollection.push_back( hit2 );
					}

				}
//				cout << "coreClusterSize : " << coreClusterSize << endl;

				if( coreClusterSize > maximumClusterSizeInCore ) {

					maximumClusterSizeInCore = coreClusterSize;
					maximumClusterSizeLayer = hit1->GetIJK().at(2);

					for( unsigned int k1=0 ; k1<corePlanCluster.size() ; k1++ ) {

						Hit *hitK1 = corePlanCluster.at(k1);
						IntVector ijk1 = hitK1->GetIJK();
						double radius = 0;
						for( unsigned int k2=0 ; k2<corePlanCluster.size() ; k2++ ) {

							Hit *hitK2 = corePlanCluster.at(k2);
							IntVector ijk2 = hitK2->GetIJK();
							double distance = (ijk1.at(0)-ijk2.at(0))* (ijk1.at(0) - ijk2.at(0)) + (ijk1.at(1)-ijk2.at(1))* (ijk1.at(1) - ijk2.at(1));

							if( sqrt(distance) > maxRadius )
								maxRadius = sqrt(distance);
						}
					}
				}
			}
//			cout << "maximumClusterSizeInCore : " << maximumClusterSizeInCore << endl;
//			cout << "maxRadius : " << maxRadius << endl;
//			cout << "maximumClusterSizeLayer : " << maximumClusterSizeLayer << endl;
//
//			cout << "firstHitLayer : " << firstHitLayer << endl;


			double coneRadius = maxRadius;

			int lastLayer = hitCol->GetLastLayer();

			double iStart = 0;
			double jStart = 0;
			int nbOfHitForStartingPoint = 0;
			for( unsigned int i=0 ; i<hitCol->size() ; i++ ) {
				if( hitCol->at(i)->GetIJK().at(2) == firstHitLayer ) {
					nbOfHitForStartingPoint++;
					iStart += hitCol->at(i)->GetIJK().at(0);
					jStart += hitCol->at(i)->GetIJK().at(1);
				}
			}
			iStart /= nbOfHitForStartingPoint;
			jStart /= nbOfHitForStartingPoint;

			ThreeVector coreStartingPoint( iStart , jStart , firstHitLayer );

			ThreeVector coneDirected( coreCog - coreStartingPoint );
			if( coneDirected.z() < 0 ) coneDirected *= -1;



			ThreeVector coneCoreStartingPoint = coreStartingPoint;

			coneCoreStartingPoint.setZ( coneCoreStartingPoint.z() -1 );

			double currentConeLength = abs(coneCoreStartingPoint.z() - maximumClusterSizeLayer);
			coneDirected.setMag( currentConeLength );
//			cout << "currentConeLength : " << currentConeLength << endl;

			Cone *bestCone = new Cone( coneCoreStartingPoint , coneDirected , coneRadius/2.0 );

//			cout << "bestCone->GetPeakPosition()" << bestCone->GetPeakPosition() << endl;

			/******************************/


			ThreeVector forwardMinimumConePeakPosition = coneCoreStartingPoint - coneBackwardDistance*coneDirected.unit();
//			cout << "forwardMinimumConePeakPosition : " << forwardMinimumConePeakPosition << endl;

			ThreeVector previousPosition = bestCone->GetPeakPosition();

//			cout << "cone peak position before looping : " << bestCone->GetPeakPosition() << endl;

			for( double coneDistance=0 ; coneDistance<coneBackwardDistance ; coneDistance+=coneStepDistance ) {


				ThreeVector currentConePeakPosition = previousPosition - coneStepDistance*coneDirected.unit();
//				cout << "coneDirected.unit() : " << coneDirected.unit() << endl;
//				cout << "currentConePeakPosition : " << currentConePeakPosition << endl;

				bestCone->SetPeakPosition( currentConePeakPosition );

				vector<IntVector> existingPadsInCone;
				vector<IntVector> touchedPadsInCone;

				this->FindPadsInCone( bestCone , existingPadsInCone , touchedPadsInCone );

				double percent = double(touchedPadsInCone.size()) / double(existingPadsInCone.size());

				if( percent < conePadsLowerLimitCut ) {

					bestCone->SetPeakPosition( previousPosition );
					break;
				}
				previousPosition = currentConePeakPosition;
			}

//			cout << "cone peak position after looping : " << bestCone->GetPeakPosition() << endl;

			core->SetStartingCone( bestCone );

			ShowerManager *showerManager = ShowerManager::GetInstance();
			Shower *shower = new Shower();
			shower->SetStartingCone( bestCone );
			shower->SetStartingPoint( bestCone->GetPeakPosition() );
			shower->AddCore( core );
			showerManager->AddShower( shower );

			TGeoManager *geoManager = SDHCALPrototype::GetGeoManager();

			if( geoManager !=0 ) {

//				ThreeVector ex(1,0,0);
//				ThreeVector ey(0,1,0);
//				ThreeVector ez(0,0,1);
//
//				double psi = bestCone->GetDirectedVector().angle( ex );
//				double phi = bestCone->GetDirectedVector().angle( ey );
//				double theta = bestCone->GetDirectedVector().angle( ez );

				ss << "cone_" << c+1;
				TGeoMaterial *vacMat = new TGeoMaterial("Vacuum",0,0,0);
				TGeoMedium *vacMed = new TGeoMedium("Vacuum",1,vacMat);
				double coneRadius = bestCone->GetRadius();
				TGeoVolume *cone = geoManager->MakeCone( ss.str().c_str() , vacMed , (currentConeLength*2.613)/2.0 ,0,0,coneRadius-coneRadius*(0.01/100.0),coneRadius);
				TGeoCone *coneShape = ( TGeoCone* )cone->GetShape();
				ThreeVector coneTrans = bestCone->GetPeakPosition() + currentConeLength/2.0 * bestCone->GetDirectedVector().unit();
				TGeoTranslation *translation = new TGeoTranslation( coneTrans.x() , coneTrans.y() , 2.613*coneTrans.z() );
				geoManager->GetTopVolume()->AddNode( cone , c+1 , translation );
				ss << "_" << c+1;
				TGeoNode *n = geoManager->GetTopVolume()->GetNode(ss.str().c_str());
				nodeList.push_back(n);
				ss.str("");
			}

		}

	}


	void ConeBeginningAlgorithm::End() {

	}


	void ConeBeginningAlgorithm::FindPadsInCone( Cone *cone , vector<IntVector> &existingPadsInCone, vector<IntVector> &touchedPadsInCone ) {

		HitManager* hitManager = HitManager::GetInstance();
		double r = cone->GetRadius();
		double l = r / sin(cone->GetTheta());

		int xMin = int(cone->GetPeakPosition().x()) - l - 1;
		int xMax = int(cone->GetPeakPosition().x()) + l + 1;
		int yMin = int(cone->GetPeakPosition().y()) - l - 1;
		int yMax = int(cone->GetPeakPosition().y()) + l + 1;
		int zMin = int(cone->GetPeakPosition().z()) - 1;
		int zMax = zMin + int(l) + 2;

//		cout << "xMin : " << xMin << endl;
//		cout << "xMax : " << xMax << endl;
//		cout << "yMin : " << yMin << endl;
//		cout << "yMax : " << yMax << endl;
//		cout << "zMin : " << zMin << endl;
//		cout << "zMax : " << zMax << endl;

		for( int i=xMin ; i<=xMax ; i++ ) {
			for( int j=yMin ; j<=yMax ; j++ ) {
				for( int k=zMin ; k<=zMax ; k++ ) {

					if( !hitManager->PadExists(i,j,k) )
						continue;
					IntVector ijk;
					ThreeVector ijkVec;
					ijkVec.set(i,j,k);
					ijk.push_back(i); ijk.push_back(j); ijk.push_back(k);

					if( cone->Contains( ijkVec ) )
						existingPadsInCone.push_back( ijk );
					else continue;

					if( !hitManager->PadIsTouched(i,j,k) )
						continue;
					touchedPadsInCone.push_back( ijk );
				}
			}
		}
	}

	int ConeBeginningAlgorithm::NbOfCoreHitsInCone( Cone *cone , Core *core ) {

		int nbOfCoreHitsInCone = 0;
		HitCollection *coreHitCollection = core->GetHitCollection();

		for( unsigned int i=0 ; i<coreHitCollection->size() ; i++ ) {

			IntVector ijk = coreHitCollection->at(i)->GetIJK();
			ThreeVector ijkVec( ijk.at(0) , ijk.at(1) , ijk.at(2) );
			if( cone->Contains( ijkVec ) )
				nbOfCoreHitsInCone++;
		}

		return nbOfCoreHitsInCone;
	}



}  // namespace 

