/*
 *
 * //DrawEvent.cc event main file
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of this software and preparation of derivative works
 * based upon this software are permitted. Any copy of this software or
 * of any derivative work must include this paragraph.
 *
 * Written by : R. Eté
 *
 */


// std includes
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <vector>
#include <stdexcept>

// baboon includes
#include "Config/SdhcalConfig.hh"
#include "Managers/HitManager.hh"
#include "Objects/Hit.hh"
#include "Utilities/Globals.hh"
#include "MarlinProcessor/CaloHitCreator.hh"

// baboon gui includes
#include "Gui/BaboonMainFrame.hh"

// lcio includes
#include "IO/LCReader.h"
#include "Exceptions.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCCollection.h"
#include "EVENT/CalorimeterHit.h"
#include "UTIL/CellIDDecoder.h"

// root includes
#include "TRint.h"
#include "TEveManager.h"
#include "TEveBox.h"
#include "TEveBrowser.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TStyle.h"
#include "TView3D.h"
#include "TRootBrowser.h"


// tclap includes
#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraint.h"

using namespace std;
using namespace baboon;
//using namespace cfgparser;

// for debug messages (devel)
#define __DEBUG__ 0

extern TEveManager *gEve;

int main (int argc ,char *argv[]) {

	/********************************************
	 * Grab the most useful environment variable
	 ********************************************/

	string pathToBaboon = getenv("BABOON_HOME");
	if( pathToBaboon.empty() )
		throw runtime_error("'BABOON_HOME' env variable is not set.\n "
		"Please source init_Baboon.sh before running.");

	SdhcalConfig *config = SdhcalConfig::GetInstance();
	config->LoadFile( pathToBaboon + string("/config/SDHCAL.cfg") );


	/*********************************
	 * Define the command line parser
	 *********************************/

	string footer = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine cmd(footer, ' ', "1.0");
	TCLAP::ValueArg<std::string> slcioFileArg(
				"f"
				,"slcio-file"
				,"slcio sdhcal event file to display"
				,false
				,""
				,"string" );
	TCLAP::ValueArg<std::string> gearFileArg(
				"g"
				,"gear-file"
				,"gear file for geometry"
				,false
				,""
				,"string" );

	cmd.add( slcioFileArg );
	cmd.add( gearFileArg );
	cmd.parse( argc, argv );


    TRint *baboonApp = new TRint("Baboon event display",0,0);

    TEveManager::Create();
    gEve->GetMainWindow()->SetWindowName("Baboon Event Display");


    /*
     * Start embedding
     */
	gEve->GetBrowser()->StartEmbedding(TRootBrowser::kLeft);

    BaboonMainFrame *mainFrame = new BaboonMainFrame();

    mainFrame->MapSubwindows();
    mainFrame->Resize();
    mainFrame->MapWindow();

    gEve->GetBrowser()->StopEmbedding("Baboon options");

    /*
     * Stop embedding
     */

    string slcioFile;
    IO::LCReader *lcReader = LCFactory::getInstance()->createLCReader() ;

    if( !slcioFileArg.getValue().empty() ) {

    	slcioFile = slcioFileArg.getValue();
    	try{
    		lcReader->open( slcioFile );
    	} catch (IO::IOException& e ) {
    		cerr << e.what() << endl;
    		cerr << "BaboonGui have to stop. Sorry!" << endl;
    		exit(1);
    	}

    	EVENT::LCCollection *lcCollection = 0;
    	EVENT::LCEvent *event = 0;
		try {
			event = lcReader->readNextEvent();
//			lcCollection = lcReader->readNextEvent()->getCollection("HCALBarrel");
		} catch (EVENT::DataNotAvailableException& e) {
			cerr << e.what() << endl;
			cerr << "BaboonGui have to stop. Sorry!" << endl;
			exit(1);
		}

		map<Hit*,TEveBox*> calorimeterHitBoxes;

		HitManager *hitManager = HitManager::GetInstance();

		CaloHitCreator *creator = new CaloHitCreator();
//		caloHitCreator->SetDecoderString( "M:3,S-1:3,I:9,J:9,K-1:6" );
		creator->SetCollectionName( "HCALBarrel" );
		creator->CreateCaloHits( event );
//		hitManager->BeginOfEvent( lcCollection );
//		hitManager->BeginOfEvent( lcCollection );
		HitCollection *hitCollection = hitManager->GetHitCollection();

		for( unsigned int i=0 ; i<hitCollection->size() ; i++ ) {

			Hit *hit = hitCollection->at(i);
			TEveBox* calorimeterHitBox = new TEveBox;
			calorimeterHitBoxes[ hit ] = calorimeterHitBox;

			if( hit->GetThreshold() == fThreshold3 )
				calorimeterHitBox->SetMainColor( kRed );
			else if( hit->GetThreshold() == fThreshold1 )
				calorimeterHitBox->SetMainColor( kBlue );
			else if( hit->GetThreshold() == fThreshold2 )
				calorimeterHitBox->SetMainColor( kGreen );
			calorimeterHitBox->SetMainTransparency( 0 );

			IntVector ijk = hit->GetIJK();
			double boxSizeFactor = 2.0 / 8.0;

			calorimeterHitBox->SetVertex(0, ijk.at(0) - boxSizeFactor , ijk.at(1) - boxSizeFactor , ijk.at(2) - boxSizeFactor );
			calorimeterHitBox->SetVertex(1, ijk.at(0) + boxSizeFactor , ijk.at(1) - boxSizeFactor , ijk.at(2) - boxSizeFactor );
			calorimeterHitBox->SetVertex(2, ijk.at(0) + boxSizeFactor , ijk.at(1) + boxSizeFactor , ijk.at(2) - boxSizeFactor );
			calorimeterHitBox->SetVertex(3, ijk.at(0) - boxSizeFactor , ijk.at(1) + boxSizeFactor , ijk.at(2) - boxSizeFactor );
			calorimeterHitBox->SetVertex(4, ijk.at(0) - boxSizeFactor , ijk.at(1) - boxSizeFactor , ijk.at(2) + boxSizeFactor );
			calorimeterHitBox->SetVertex(5, ijk.at(0) + boxSizeFactor , ijk.at(1) - boxSizeFactor , ijk.at(2) + boxSizeFactor );
			calorimeterHitBox->SetVertex(6, ijk.at(0) + boxSizeFactor , ijk.at(1) + boxSizeFactor , ijk.at(2) + boxSizeFactor );
			calorimeterHitBox->SetVertex(7, ijk.at(0) - boxSizeFactor , ijk.at(1) + boxSizeFactor , ijk.at(2) + boxSizeFactor );

			gEve->AddElement( calorimeterHitBox );

		}

		gEve->Redraw3D(kTRUE);


    }

    gEve->FullRedraw3D(kTRUE);

	SdhcalConfig::Kill();
	baboonApp->Run();

//	delete baboonApp;

	return 0;
}



