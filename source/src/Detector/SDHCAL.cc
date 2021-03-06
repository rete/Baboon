  /// \file SDHCAL.cc
/*
 *
 * SDHCAL.cc source template generated by fclass
 * Creation date : ven. juil. 12 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Detector/SDHCAL.hh"

using namespace std;

namespace baboon {

	SDHCAL::SDHCAL( const std::string &instanceName )
		: Calorimeter( instanceName , kHcalEndcap ) {

		viewContent = true;
	}

	SDHCAL::~SDHCAL() {

	}


	Return SDHCAL::ReadSettings( const gear::CalorimeterParameters *caloParameters ) {

		try {

			cellSize0 = caloParameters->getLayerLayout().getCellSize0(0);
			cellSize1 = caloParameters->getLayerLayout().getCellSize1(0);
			absorberThickness = caloParameters->getLayerLayout().getAbsorberThickness(0);
			layerThickness = caloParameters->getLayerLayout().getThickness(0);
			numberOfLayers = caloParameters->getLayerLayout().getNLayers();
			repeatX = caloParameters->getIntVal("repeatX");
			repeatY = caloParameters->getIntVal("repeatY");

		}
		catch ( std::exception &e ) {
			std::cerr << "Couln't load SDHCAL detector gear parameters!" << std::endl;
			std::cerr << e.what() << std::endl;
			return BABOON_ERROR("Couln't read correctly the settings");
		}

		// Get the radiation length.
		try {
			radiationLength = caloParameters->getDoubleVal("radiationLength");
		}
		catch ( std::exception &e ) {
			std::cout << "Radiation length not found!" << std::endl;
			radiationLength = 0;
		}

		// Get the interaction length.
		try {
			interactionLength = caloParameters->getDoubleVal("interactionLength");
		}
		catch ( std::exception &e ) {
			std::cout << "Interaction length not found!" << std::endl;
			interactionLength = 0;
		}

		// Get the global X shift length.
		try {
			globalShiftX = caloParameters->getDoubleVal("globalShiftX");
		}
		catch ( std::exception &e ) {
			std::cout << "Global X shift not found!" << std::endl;
			globalShiftX = 0.0;
		}

		// Get the global Y shift length.
		try {
			globalShiftY = caloParameters->getDoubleVal("globalShiftY");
		}
		catch ( std::exception &e ) {
			std::cout << "Global Y shift not found!" << std::endl;
			globalShiftY = 0.0;
		}

		// Get the global Z shift length.
		try {
			globalShiftZ = caloParameters->getDoubleVal("globalShiftZ");
		}
		catch ( std::exception &e ) {
			std::cout << "Global Z shift not found!" << std::endl;
			globalShiftZ = 0.0;
		}

		return BABOON_SUCCESS();
	}


	Return SDHCAL::BuildGeometry( TGeoManager *manager , TGeoVolume *topVolume ) {

		if( this->isBuilt )
			return BABOON_ALREADY_INITIALIZED("SDHCAL already built!");


		for( unsigned int l=0 ; l<numberOfLayers ; l++ ) {

			TGeoMaterial *vacuumMaterial = new TGeoMaterial("Vacuum", 0, 0, 0); // dummy material
			TGeoMedium *vacuum = new TGeoMedium("Vacuum",1, vacuumMaterial);

			std::ostringstream oss;
			oss << "Layer_" << l ;
			TGeoVolume *layer = manager->MakeBox(oss.str().c_str(),vacuum,repeatX*cellSize0/2.0,repeatY*cellSize1/2.0,absorberThickness/2.0);

			layer->SetLineColor(kGray);
			layer->SetTransparency(95);
			layer->SetVisibility(true);

			topVolume->AddNode(layer, 0, new TGeoTranslation(0, 0, absorberThickness/2.0 + l*layerThickness ) );
		}

		isBuilt = true;
		return BABOON_SUCCESS();
	}



	Return SDHCAL::ViewCalorimeterContent( TEveElement *parent , const HitDisplayMode displayMode ) {

		TEveBoxSet *hits = new TEveBoxSet("SDHCAL_HITS");
		hits->Reset(TEveBoxSet::kBT_AABox, kTRUE, 64);
		hits->SetOwnIds(kTRUE);
		hits->SetPickable(kTRUE);
		hits->SetMainTransparency(100);
		hits->SetAntiFlick(kTRUE);

		for( unsigned int i=0 ; i<caloHitCollection->size() ; i++ ) {

			CaloHit *caloHit = caloHitCollection->at(i);
			IntVector ijk = caloHit->GetIJK();

			double posX , posY , posZ = 0.;

			posX = ijk.at(0)*cellSize0 - repeatX*cellSize0/2.0 - cellSize0/2.0;
			posY = ijk.at(1)*cellSize1 - repeatY*cellSize1/2.0 - cellSize1/2.0;
			posZ = ijk.at(2)*layerThickness + absorberThickness + (layerThickness-absorberThickness)/2.0;

			hits->AddBox( posX , posY , posZ , cellSize0*0.8 , cellSize1*0.8 , layerThickness-absorberThickness );
			int color = BaboonMonitoring::GetInstance()->GetCaloHitColor( caloHit , displayMode );

			// transparent hit ( not display )
			if( color < 0 )
				hits->DigitColor( color , 100 );
			else
				hits->DigitColor( color , 0 );
		}

		if( parent ) {
			parent->AddElement( hits );
		}
		else {
			gEve->AddElement( hits );
			gEve->Redraw3D();
		}

		return BABOON_SUCCESS();
	}

}  // namespace 

