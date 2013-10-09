  /// \file Calorimeter.cc
/*
 *
 * Calorimeter.cc source template generated by fclass
 * Creation date : mer. juin 12 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Detector/Calorimeter.hh"

namespace baboon {

	Calorimeter::Calorimeter( const std::string &detectorName , const CalorimeterType &caloType )
		: Detector( detectorName , kCalorimeter ),
		  calorimeterType(caloType),
		  cellSize0(0),
		  cellSize1(0),
		  repeatX(0),
		  repeatY(0),
		  absorberThickness(0),
		  layerThickness(0),
		  numberOfLayers(0),
		  globalShiftX(0),
		  globalShiftY(0),
		  globalShiftZ(0),
		  interactionLength(0),
		  radiationLength(0),
		  hitDisplayMode( kDisplayThresholds ) {

		caloHitCollection = new CaloHitCollection();
	}



	Calorimeter::~Calorimeter() {

		for( unsigned int i=0 ; caloHitCollection->size() ; i++ )
			if( caloHitCollection->at(i) != 0 )
				delete caloHitCollection->at(i);

		caloHitCollection->clear();
		delete caloHitCollection;
	}



	Return Calorimeter::AddCaloHit( CaloHit *caloHit ) {

		if( !this->isBuilt )
			return BABOON_ERROR("Calorimeter is not initialized. Can't add a hit");

		if( caloHitCollection->empty() ) {

			BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , caloHitCellHelper.AddCaloHit( caloHit ) );
			caloHitCollection->push_back( caloHit );
			return BABOON_SUCCESS();
		}
		CaloHitCollection::iterator caloHitIt = std::find( caloHitCollection->begin() , caloHitCollection->end() , caloHit );

		if( caloHitIt == caloHitCollection->end() ) {
			caloHitCollection->push_back( caloHit );
			BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , caloHitCellHelper.AddCaloHit( caloHit ) );
			return BABOON_SUCCESS();
		}
		return BABOON_ALREADY_PRESENT("Calo hit already registered in the calorimeter");
	}



	Return Calorimeter::RemoveCaloHit( CaloHit *caloHit ) {

		if( !this->isBuilt )
			return BABOON_ERROR("Calorimeter is not initialized. Can't add a hit");

		if( caloHitCollection->empty() )
			return BABOON_INVALID_PARAMETER("Calo hit collection is empty!");

		CaloHitCollection::iterator caloHitIt = std::find( caloHitCollection->begin() , caloHitCollection->end() , caloHit );

		if( caloHitIt != caloHitCollection->end() ) {

			caloHitCollection->erase( caloHitIt );
			BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , caloHitCellHelper.RemoveCaloHit( caloHit ) );
			return BABOON_SUCCESS("Calo hit correctly removed");
		}

		return BABOON_NOT_FOUND("Calo hit was not registered in the calorimeter");
	}



	Return Calorimeter::ReadSettings( const gear::GearParameters *parameters ) {

		try {
			const gear::CalorimeterParameters *caloParams = dynamic_cast<const gear::CalorimeterParameters *>(parameters);
			BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , this->ReadSettings( caloParams ) );
		}
		catch( std::exception &e ) {
			return BABOON_ERROR("While reading settings : wrong parameters type! e.what() : " + std::string(e.what()));
		}
		return BABOON_SUCCESS();
	}


	Return Calorimeter::ClearContent() {

		caloHitCellHelper.Clear();

		if( caloHitCollection != 0 ) {
			for( unsigned int i=0 ; i<caloHitCollection->size() ; i++ ) {

				if( caloHitCollection->at(i) != 0 ) {
					delete caloHitCollection->at(i);
					caloHitCollection->at(i) = 0;
				}
				else
					return BABOON_ERROR("Hit is a NULL pointer!");
			}
			caloHitCollection->clear();
		}
		else
			return BABOON_ERROR("Hit collection is a NULL pointer!");

		return BABOON_SUCCESS("Content cleared!");
	}



	CaloHit *Calorimeter::GetCaloHitAt( unsigned int I , unsigned int J , unsigned int K ) {

		if( this->IsPadFired( I , J , K ) )
			return caloHitCellHelper.GetCaloHitAt( I , J , K );

		else return 0;
	}



	bool Calorimeter::IsPadFired( unsigned int I , unsigned int J , unsigned int K ) {

		return caloHitCellHelper.IsPadFired( I , J , K );
	}


}  // namespace 

