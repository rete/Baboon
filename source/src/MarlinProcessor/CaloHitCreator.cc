  /// \file CaloHitCreator.cc
/*
 *
 * CaloHitCreator.cc source template generated by fclass
 * Creation date : jeu. mai 9 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "MarlinProcessor/CaloHitCreator.hh"


using namespace std;
using namespace EVENT;

namespace baboon {

	CaloHitCreator::CaloHitCreator() {

		collectionName = "HCALBarrel";
		decoderString = "M:3,S-1:3,I:9,J:9,K-1:6";
	}

	CaloHitCreator::~CaloHitCreator() {
	}

	Return CaloHitCreator::CreateCaloHits( LCEvent *evt ) {


		try {

			LCCollection *collection = evt->getCollection( collectionName );
			UTIL::CellIDDecoder<CalorimeterHit>::setDefaultEncoding( decoderString );
			UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder( collection );

			for( unsigned int i=0 ; i<collection->getNumberOfElements() ; i++ ) {

				CalorimeterHit *caloHit = static_cast<CalorimeterHit*> ( collection->getElementAt(i) );
				HitParameters hitParams;
				ThreeVector position;
				position.setX( caloHit->getPosition()[0] );
				position.setY( caloHit->getPosition()[1] );
				position.setZ( caloHit->getPosition()[2] );
				hitParams.position = position;
				try {
	//				cout << "hit energy : " << caloHit->getEnergy() << endl;
					if( caloHit->getEnergy() == 1.0 ) hitParams.threshold = fThreshold1;
					else if( caloHit->getEnergy() == 2.0 ) hitParams.threshold = fThreshold2;
					else if( caloHit->getEnergy() == 3.0 ) hitParams.threshold = fThreshold3;
					else throw ThresholdUndefinedException();
				}
				catch ( ThresholdUndefinedException &e ) {

					stringstream ss;
					ss << "Undefined threshold exception thrown. Energy is not unit of SDHCAL threshold (1,2,3)" << endl;
					cerr << ss.str();
					return S_ERROR( ss.str() );
				}
				int I = cellIdDecoder(caloHit)["I"];
				int J = cellIdDecoder(caloHit)["J"];
				int K = cellIdDecoder(caloHit)["K-1"];
				IntVector ijk;
				ijk.push_back( I );
				ijk.push_back( J );
				ijk.push_back( K );
				hitParams.ijk = ijk;

				HitManager::GetInstance()->RegisterNewHit( hitParams );
			}
		}
		catch ( DataNotAvailableException &e ) {

			stringstream ss;
			ss << "LCIO exception thrown : " << e.what() << endl;
			cerr << ss.str();
			return S_ERROR( ss.str() );
		}
		catch ( Exception &e ) {

			stringstream ss;
			ss << "Failed to create hit. " << e.what() << endl;
			cerr << ss.str();
			return S_ERROR( ss.str() );
		}

		return S_OK();
	}



}  // namespace



