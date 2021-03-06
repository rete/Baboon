  /// \file BaboonMonitoring.cc
/*
 *
 * BaboonMonitoring.cc source template generated by fclass
 * Creation date : jeu. juil. 11 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Monitoring/BaboonMonitoring.hh"

#include "TCanvas.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TTree.h"

#include <TEveManager.h>
#include <TEveEventManager.h>
#include <TEveViewer.h>
#include <TEvePointSet.h>
#include <TEveArrow.h>
#include <TEveRGBAPalette.h>
#include <TGLViewer.h>
#include <TEveBrowser.h>

#include <TGeoXtru.h>
#include <TEveGeoShapeExtract.h>
#include <TEveGeoShape.h>

#include <TGeometry.h>
#include <TGeoMaterial.h>
#include <TGeoManager.h>
#include <TEveGeoNode.h>

#include "TParticlePDG.h"

#include <TSystem.h>
#include <TGeoXtru.h>
#include <TGeoMatrix.h>
#include <TGeoTube.h>
#include <TEveBoxSet.h>
#include <TGeoCompositeShape.h>
#include <TEveScene.h>
#include <TGLLightSet.h>

#include <vector>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <cmath>
#include <fcntl.h>
#include <limits>
#include <algorithm>


using namespace std;


namespace baboon {


	BaboonMonitoring *BaboonMonitoring::instance = 0;
	bool BaboonMonitoring::openEveEvent = false;
	int  BaboonMonitoring::eventDisplayCounter = 0;
	bool BaboonMonitoring::isTEveInitialized = false;
	bool BaboonMonitoring::isEnable = false;
	std::map<string,int> BaboonMonitoring::tagToColorMap;



	BaboonMonitoring::BaboonMonitoring() {

		rootApplication = new TApplication("BaboonMonitoring", 0, 0);
	    rootApplication->SetReturnFromRun(kTRUE);

	    eveManager = TEveManager::Create( false );
	}

	BaboonMonitoring::~BaboonMonitoring() {

		tagToColorMap.clear();

		if( isTEveInitialized ) {

			TEveManager::Terminate();
			gSystem->ProcessEvents();
			eveManager = 0;
		}
	    rootApplication->Terminate(0);
	}


	BaboonMonitoring *BaboonMonitoring::GetInstance() {

		if( instance == 0 ) {

			instance = new BaboonMonitoring();
	        TColor::CreateColorWheel();

	        // Build the tag to color map
	        tagToColorMap[ UndefinedTag::Tag() ] = kWhite;
	        tagToColorMap[ TrackTag::Tag() ] = kRed;
	        tagToColorMap[ TrackExtremityTag::Tag() ] = kGreen;
	        tagToColorMap[ IsolatedTag::Tag() ] = kBlue;
	        tagToColorMap[ CoreTag::Tag() ] = kMagenta;
	        tagToColorMap[ NoiseTag::Tag() ] = kYellow;

	        gStyle->SetPalette(1);
	        gStyle->SetNumberContours(99);
		}
		return instance;
	}


	void BaboonMonitoring::Kill() {

		delete GetInstance();
	}

	void BaboonMonitoring::AddElement( TEveElement *element , TEveElement *parent ) {

		InitializeEve();
		if( eveManager && isEnable ) {
			eveManager->AddElement( element , parent );
		}
	}

	void BaboonMonitoring::Pause() const {

		if( !isEnable )
			return;

	#ifdef __unix__
	    std::cout << "Press return to continue ..." << std::endl;
	    int flag = fcntl(1, F_GETFL, 0);

	    int key = 0;
	    while(true)
	    {
	        gSystem->ProcessEvents();
	        fcntl(1, F_SETFL, flag | O_NONBLOCK);
	        key = getchar();

	        if((key == '\n') || (key == '\r'))
	            break;

	        usleep(1000);
	    }

	    fcntl(1, F_SETFL, flag);
	#else
	    std::cout << "BaboonMonitoring::Pause() is only implemented for unix operating systems." << std::endl;
	#endif
	}



	void BaboonMonitoring::InitializeEve( Char_t transparency ) {

		if( !isEnable )
			return;

	    std::stringstream sstr;
	    sstr << "Event Display " << eventDisplayCounter;

	    // If Eve instantiated
	    if ( isTEveInitialized )
	    {
	        TEveEventManager *currentEvent = gEve->GetCurrentEvent();

	        if (currentEvent)
	            currentEvent->SetElementNameTitle( sstr.str().c_str() , sstr.str().c_str() );


	        if (!openEveEvent)
	        {
	        	cout << "Event created" << endl;
	            gEve->AddEvent( new TEveEventManager( sstr.str().c_str(),sstr.str().c_str() ) );
	            openEveEvent = true;
	            eventDisplayCounter++;
	        }

	        return;
	    }


	    // If Eve is not instantiated
//	    gSystem->Load("libGeom");

	    DetectorManager *man = DetectorManager::GetInstance();

	    if( !man->geometryBuilt )
	    	BABOON_THROW_RESULT( BABOON_NOT_INITIALIZED("Detector geometry is not built! Couldn't display event!") );


	    try
	    {
	        std::cout << "BaboonMonitoring::InitializeEve(): ";
	        const char *display(::getenv("DISPLAY"));

	        if (NULL == display)
	        {
	            std::cout << "DISPLAY environment not set" << std::endl;
	        }
	        else
	        {
	            std::cout << "DISPLAY environment set to " << display << std::endl;
	        }

	        TEveManager::Create();
	    }
	    catch (TEveException &tEveException)
	    {
	        std::cout << "BaboonMonitoring::InitializeEve(): Caught TEveException: " << tEveException.what() << std::endl;

	        try
	        {
	            std::cout << "BaboonMonitoring::InitializeEve(): Attempt to release ROOT from batch mode." << std::endl;
	            gROOT->SetBatch(kFALSE);
	            TEveManager::Create();
	        }
	        catch (TEveException &tEveException)
	        {
	            std::cout << "BaboonMonitoring::InitializeEve(): Caught TEveException: " << tEveException.what() << std::endl;
	            throw std::exception();
	        }
	    }
	    catch ( std::exception &e ) {
	    	std::cerr << "BaboonMonitoring::InitializeEve(): Caught std::exception: " << e.what() << std::endl;
	    	throw std::exception();
	    }

	    TGeoNode *geoNode = man->geoManager->GetTopNode();
	    TEveGeoTopNode *eveGeoTopNode = new TEveGeoTopNode( man->geoManager, geoNode );
	    eveGeoTopNode->SetVisLevel(1);
	    eveGeoTopNode->GetNode()->GetVolume()->SetVisibility(kFALSE);

	    gEve->AddGlobalElement( eveGeoTopNode );
	    gEve->GetBrowser()->MapWindow();

	    TGLViewer *viewerGL = gEve->GetDefaultGLViewer();
	    viewerGL->ColorSet().Background().SetColor(kBlack);
	    viewerGL->GetLightSet()->SetUseSpecular(kFALSE);

	    gEve->GetBrowser()->GetMainFrame()->SetWindowName("Baboon Monitoring");
	    gEve->AddEvent(new TEveEventManager(sstr.str().c_str(),sstr.str().c_str()));
	    eventDisplayCounter++;

	    isTEveInitialized = true;
	    openEveEvent = true;
	}


	int BaboonMonitoring::GetColorBetweenMinAndMax( double min , double max , double value ) {

		return (51 + ( (value - min)*( 100 - 51 ) )/( max - min ) );
	}

	//------------------------------------------------------------------------------------------------------------------------------------------



	void BaboonMonitoring::ViewDetectorsContent() {

		DetectorManager *man = DetectorManager::GetInstance();

		StringVector detectorList = man->GetDetectorList();

		for( unsigned int d=0 ; d<detectorList.size() ; d++ ) {

			Detector *det = man->GetDetector( detectorList.at(d) );
			if( det->GetDetectorType() == kCalorimeter ) {
				Calorimeter *calo = (Calorimeter *)det;
				if( calo->IsViewContentSet() ) {
					HitDisplayMode hitDisplayMode = calo->GetHitDisplayMode();
					calo->ViewCalorimeterContent( 0 , hitDisplayMode );
				}
			}
		}
	}


	void BaboonMonitoring::ViewEvent() {

		if( !isEnable )
			return;

	    InitializeEve();

	    this->ViewDetectorsContent();

	    gEve->Redraw3D( true );
	    this->Pause();

	    TEveEventManager* current = gEve->GetCurrentEvent();
	    if (current)
	        current->SetRnrSelfChildren(kFALSE,kFALSE);

	    openEveEvent = false;
	    std::cout << "View done" << std::endl;
	}


	int BaboonMonitoring::GetCaloHitColor( CaloHit *caloHit , const HitDisplayMode displayMode ) {

		if( displayMode == kDisplayUniform ) {
			return kOrange+2;
		}
		else if( displayMode == kDisplayThresholds ) {

			CaloHitThreshold thr = caloHit->GetThreshold();
			if( thr == fCaloHitThr1 ) {
				return kGreen;
			}
			else if( thr == fCaloHitThr2 ) {
				return kBlue;
			}
			else if( thr == fCaloHitThr3 ) {
				return kRed;
			}
			else return kWhite;
		}
		else if( displayMode == kDisplayTags ) {

			return caloHit->GetTag().GetColor();
		}
		else if( displayMode ==  kDisplayShowers ) {

			ShowerCollection *showers = ShowerManager::GetInstance()->GetShowerCollection();

			if( showers->empty() )
				return kWhite;

			for( unsigned int sh=0 ; sh<showers->size() ; sh++ ) {

				if( showers->at(sh)->Contains( caloHit ) ) {
					return sh+2;
				}
			}

			// if not found return kWhite
			return kWhite;
		}
		// kCustom case
		else {
			return caloHit->GetColor();
		}

	}


	void BaboonMonitoring::DrawConnector( const ThreeVector &position1 , const ThreeVector &position2 , int color ) {

		if( isEnable ) {

			TEveArrow *connectionArrow = new TEveArrow(
					 position2.x() - position1.x()
					,position2.y() - position1.y()
					,position2.z() - position1.z()
					,position1.x()
					,position1.y()
					,position1.z() );

			connectionArrow->SetMainColor( color );
			connectionArrow->SetPickable( true );
//			connectionArrow->SetTubeR( 0.05 );
//			connectionArrow->SetConeR( 0.008 );
//			connectionArrow->SetConeL( 0.008 );
			instance->AddElement( connectionArrow );
		}
	}


}  // namespace 

