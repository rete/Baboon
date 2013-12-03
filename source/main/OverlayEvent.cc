/*
 *
 * //OverlayEvent.cc main file
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
#include <iomanip>

// CfgParser includes
#include "CfgParser/CfgParser.hh"
#include "CfgParser/Data.hh"

// lcio includes
#include "IO/LCReader.h"
#include "IO/LCWriter.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCCollection.h"
#include "EVENT/CalorimeterHit.h"
#include "IMPL/CalorimeterHitImpl.h"
#include "UTIL/LCTOOLS.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCTOOLS.h"
#include "UTIL/CellIDDecoder.h"
#include "UTIL/CellIDEncoder.h"
#include "lcio.h"

// baboon includes
#include "Overlay/Overlayer.hh"
#include "Objects/CaloHit.hh"
#include "Utilities/ReturnValues.hh"
#include "Managers/AlgorithmManager.hh"
#include "Algorithm/OverlayEventAlgorithm.hh"
#include "Algorithm/Calorimetry/TrackFinderAlgorithm.hh"
#include "Algorithm/Calorimetry/ClusteringAlgorithm.hh"
#include "Utilities/Internal.hh"
#include "Utilities/ReturnValues.hh"
#include "Utilities/CaloHitHelper.hh"
#include "Detector/SDHCAL.hh"

//#include "gear/GearXML.h"
#include "gear/GearMgr.h"

// tclap includes
#include "tclap/CmdLine.h"

using namespace std;
using namespace cfgparser;
using namespace baboon;
using namespace EVENT;
using namespace UTIL;
using namespace gear;

// for debug messages (devel)
#define __DEBUG__ 0


int main (int argc ,char *argv[]) {

	/********************************************
	 * Grab the most useful environment variable
	 ********************************************/

	char *babHome = NULL;
	babHome = getenv("BABOON_HOME");
	if( babHome == NULL )
		throw runtime_error("'BABOON_HOME' env variable is not set.\n "
		"Please source init_Baboon.sh before running.");

	string baboonHome(babHome);

	/*********************************
	 * Define the command line parser
	 *********************************/

	string footer = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine cmd(footer, ' ', "1.0");
	TCLAP::ValueArg<std::string> cfgFileArg(
				"f"
				,"cfgfile"
				,"configuration file to run OverlayEvent"
				,true
				,""
				,"string" );
	TCLAP::SwitchArg discreteModeArg(
				"d"
				,"discrete"
				,"Set the discrete mode. No screen printing."
				,false);

	TCLAP::SwitchArg gridModeArg(
				"g"
				,"grid"
				,"Set the grid mode. Allow to run OverlayEvent on the grid without using cin checking ."
				,false);

	cmd.add( cfgFileArg );
	cmd.add( discreteModeArg );
	cmd.add( gridModeArg );
	cmd.parse( argc, argv );

	bool discreteMode = discreteModeArg.getValue();


	/*********************************************************************
	 * Open a CfgParser and read the global configuration (SDHCAL.cfg)
	 * and the configuration file to overlay events. The user must
	 * provide the last one while the first one contains default settings
	 * about SDHCAL. Do not modify this one!
	 *********************************************************************/

	vector<int> nbOfPadsXYZ;
	SdhcalConfig *config = SdhcalConfig::GetInstance();
	config->LoadFile( baboonHome + "/config/SDHCAL.cfg" );
	config->GetData("pads").GetValue("nbOfPadsXYZ",&nbOfPadsXYZ);


	// get the input configuration file from command line argv[1]
	string inputFile1 = "";
	string inputCollectionName1 = "";
	string inputFile2 = "";
	string inputCollectionName2 = "";
	string slcioOutputFile = "";
	string codingPattern1 = "";
	string codingPattern2 = "";
	string outputCollectionName = "";
	string outputCodingPattern = "";
	string gearFile = "";
	string algorithmConfigFile = "";
	int separationDistance = 0;
	int nbOfEventsToOverlay = 0;

	string overlayCfgFileName = "";

	overlayCfgFileName = cfgFileArg.getValue();


	CfgParser *parser = new CfgParser( overlayCfgFileName );
	parser->Read();

	parser->GetValue("general","gearFile",&gearFile);
	parser->GetValue("general","algorithmConfigFile",&algorithmConfigFile);

	parser->GetValue("input1","slcioFile",&inputFile1);
	parser->GetValue("input1","collectionName",&inputCollectionName1);
	parser->GetValue("input1","codingPattern",&codingPattern1);

	parser->GetValue("input2","slcioFile",&inputFile2);
	parser->GetValue("input2","collectionName",&inputCollectionName2);
	parser->GetValue("input2","codingPattern",&codingPattern2);

	parser->GetValue("output","nbOfEventsToOverlay",&nbOfEventsToOverlay);
//	parser->GetValue("output","separationDistance",&separationDistance);
	parser->GetValue("output","slcioFile",&slcioOutputFile);
	parser->GetValue("output","collectionName",&outputCollectionName);
	parser->GetValue("output","codingPattern",&outputCodingPattern);


	/****************************************
	 *          Calorimeter loading
	 ****************************************/

	GearXML gearParser( gearFile );
	GearMgr *gearMgr = gearParser.createGearMgr();

	DetectorManager * detectorMan = DetectorManager::GetInstance();
	detectorMan->RegisterDetector( new SDHCAL("SDHCAL1") );
	detectorMan->RegisterDetector( new SDHCAL("SDHCAL2") );
	BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , detectorMan->Init( gearMgr ) );


	/****************************************
	 *          Algorithm loading
	 ****************************************/

	AlgorithmManager *algoMan = AlgorithmManager::GetInstance();

	// Register overlay algo
	algoMan->RegisterAlgorithm( new OverlayEventAlgorithm() );

	// Register track finder algo
	algoMan->RegisterAlgorithm( new TrackFinderAlgorithm() );

	// Register clustering algo
	algoMan->RegisterAlgorithm( new ClusteringAlgorithm() );

	algoMan->SetConfigFileName( algorithmConfigFile );
	algoMan->Initialize();

	OverlayEventAlgorithm *overlayer = (OverlayEventAlgorithm *) algoMan->GetAlgorithm( "OverlayEventAlgorithm" );
	Data data = overlayer->GetData();
	data.GetValue("separationDistance", &separationDistance);

	/*****************************************
	 * LC readers (input) and writer (output)
	 *****************************************/

	IO::LCReader* lcReader1 = IOIMPL::LCFactory::getInstance()->createLCReader();
	IO::LCReader* lcReader2 = IOIMPL::LCFactory::getInstance()->createLCReader();
	IO::LCWriter* lcWriter  = IOIMPL::LCFactory::getInstance()->createLCWriter();
	lcReader1->open( inputFile1 );
	lcReader2->open( inputFile2 );
	lcWriter->open( slcioOutputFile , LCIO::WRITE_NEW );
	EVENT::LCEvent *evt1( nullptr );
	EVENT::LCEvent *evt2( nullptr );

	/*******************************************
	 * Check if the number of events to overlay
	 * is compatible with input files.
	 *******************************************/

	int nbOfEvents1 = lcReader1->getNumberOfEvents();
	int nbOfEvents2 = lcReader2->getNumberOfEvents();
	int maxNbOfEvents = min(nbOfEvents1,nbOfEvents2);
	bool gridMode = gridModeArg.getValue();


	if( nbOfEvents1 < nbOfEventsToOverlay || nbOfEvents2 < nbOfEventsToOverlay ) {
		if(!discreteMode) cout << "  WARNING!! : Inuput files don't have enough event to overlay.\n"
		"              A maximum of " << min(nbOfEvents1,nbOfEvents2) << " will be overlapped.\n";

		// cin get the response given by the user to know if
		// the program have to continue or not. Disable it with option -g ()
		if(!gridMode) {
			string response;
			while(response != "yes" && response != "no") {
				cout <<"Do you want to continue? (yes/no) : ";
				response.clear();
				cin >> response;
			}
			if(response == "no") {
				exit(0);
			}
		}
		nbOfEventsToOverlay = min(nbOfEvents1,nbOfEvents2);
	}


	/********************************************************************
	 * Loop over events. Overlap two events with two input files.
	 * If by redefining a hit position the hit is out of the sdhcal
	 * the hit is deleted from the output collection. The list of skipped
	 * event is printed at the end for the user. Stop overlaying when the
	 * a given number of events are overlaid (see cfg file).
	 ********************************************************************/

	int nbOfSkippedEvents = 0;
	int nbOfOverlaidEvents = 0;
	int evtID=0;
	vector<int> skippedEventsVec;
	int lastPercentDisplayed = 0;
	int nbOfLostHits = 0;
	nbOfEventsToOverlay = min(maxNbOfEvents,nbOfEventsToOverlay);


	while( nbOfOverlaidEvents != nbOfEventsToOverlay ) {

		// Progress display
		if(!discreteMode) {
			double percent = (nbOfOverlaidEvents*1.0) / (1.0*nbOfEventsToOverlay) * 10;
			if( round(percent) > lastPercentDisplayed ) {
				lastPercentDisplayed = round(percent);
				cout << "Processing ... " << lastPercentDisplayed*10 << "%" << endl;
			}
		}

		evt1 = lcReader1->readNextEvent();
		evt2 = lcReader2->readNextEvent();


		/*****************************************************************
		 * Make a copy of the input collections in HitCollection type
		 * Each collection is tagged via hit->setType(1) or hit->setType(2)
		 ******************************************************************/

		EVENT::LCCollection *lcCollection1 = evt1->getCollection(inputCollectionName1);
		EVENT::LCCollection *lcCollection2 = evt2->getCollection(inputCollectionName2);

		CaloHitCollection *collection1 = new CaloHitCollection();
		CaloHitCollection *collection2 = new CaloHitCollection();

		// Copy collection 1
		UTIL::CellIDDecoder<CalorimeterHit>::setDefaultEncoding( codingPattern1 );
		UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder1( lcCollection1 );

		for(int hitID=0 ; hitID<lcCollection1->getNumberOfElements() ; hitID++) {

			EVENT::CalorimeterHit *hit = static_cast<EVENT::CalorimeterHit *> (lcCollection1->getElementAt(hitID) );

			CaloHit *caloHit = new CaloHit( fSemiDigitalCaloHit );

			ThreeVector position;
			position.setX( hit->getPosition()[0] );
			position.setY( hit->getPosition()[1] );
			position.setZ( hit->getPosition()[2] );
			caloHit->SetPosition( position );

			int I = cellIdDecoder1( hit )["I"];
			int J = cellIdDecoder1( hit )["J"];
			int K = cellIdDecoder1( hit )["K-1"];

			caloHit->SetIJK( I , J , K );
			caloHit->SetTypeID(1);
			caloHit->SetTime( hit->getTime() );

			if( hit->getEnergy() == 1.0 ) caloHit->SetThreshold( fCaloHitThr1 );
			else if( hit->getEnergy() == 2.0 ) caloHit->SetThreshold( fCaloHitThr2 );
			else if( hit->getEnergy() == 3.0 ) caloHit->SetThreshold( fCaloHitThr3 );
			else throw runtime_error("Calo hit energy is not 1.0 , 2.0 or 3.0 as expected for SDHCAL. Check your inputs!");

			collection1->push_back( caloHit );
		}

		// Copy collection 2
		UTIL::CellIDDecoder<CalorimeterHit>::setDefaultEncoding( codingPattern2 );
		UTIL::CellIDDecoder<CalorimeterHit> cellIdDecoder2( lcCollection2 );

		for(int hitID=0 ; hitID<lcCollection2->getNumberOfElements() ; hitID++) {

			EVENT::CalorimeterHit *hit = static_cast<EVENT::CalorimeterHit *> (lcCollection2->getElementAt(hitID) );

			CaloHit *caloHit = new CaloHit( fSemiDigitalCaloHit );

			ThreeVector position;
			position.setX( hit->getPosition()[0] );
			position.setY( hit->getPosition()[1] );
			position.setZ( hit->getPosition()[2] );
			caloHit->SetPosition( position );

			int I = cellIdDecoder2( hit )["I"];
			int J = cellIdDecoder2( hit )["J"];
			int K = cellIdDecoder2( hit )["K-1"];

			caloHit->SetIJK( I , J , K );
			caloHit->SetTypeID(2);
			caloHit->SetTime( hit->getTime() );

			if( hit->getEnergy() == 1.0 ) caloHit->SetThreshold( fCaloHitThr1 );
			else if( hit->getEnergy() == 2.0 ) caloHit->SetThreshold( fCaloHitThr2 );
			else if( hit->getEnergy() == 3.0 ) caloHit->SetThreshold( fCaloHitThr3 );
			else throw runtime_error("Calo hit energy is not 1.0 , 2.0 or 3.0 as expected for SDHCAL. Check your inputs!");

			collection2->push_back( caloHit );
		}

		OverlayEventAlgorithm *overlayer = (OverlayEventAlgorithm *) algoMan->GetAlgorithm( "OverlayEventAlgorithm" );
		overlayer->SetCollection1( collection1 );
		overlayer->SetCollection2( collection2 );
		SDHCAL *sdhcal1 = (SDHCAL *) detectorMan->GetDetector("SDHCAL1");
		SDHCAL *sdhcal2 = (SDHCAL *) detectorMan->GetDetector("SDHCAL2");
		overlayer->SetCalorimeter1( sdhcal1 );
		overlayer->SetCalorimeter2( sdhcal2 );

		overlayer->Process();

//		// Compute the centers of gravity
//		ThreeVector cog1(0,0,0);
//		BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , CaloHitHelper::ComputeBarycenter( collection1 , cog1 ) );
//
//		ThreeVector cog2(0,0,0);
//		BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , CaloHitHelper::ComputeBarycenter( collection2 , cog2 ) );
//
//		Overlayer overlayer( collection1 , collection2 );

		// Separate the two events by the given separation distance
		// and center it at the center of the SDHCAL in the x direction
		// and re-center them on the y axis
//		overlayer.SetTranslations(
//				new ThreeVector( nbOfPadsXYZ.at(0)/2.0 + separationDistance / 2.0 - cog1.x()
//								 ,nbOfPadsXYZ.at(1)/2.0 - cog1.y()
//								,0)
//			   ,new ThreeVector( nbOfPadsXYZ.at(0)/2.0 - separationDistance / 2.0 - cog2.x()
//					   	   	    ,nbOfPadsXYZ.at(1)/2.0 - cog2.y()
//					   	   	    ,0) );
//
//		overlayer.OverlayCollections();



		// Grab the final overlaid collection and fill more
		// information in the collection and in the event.
//		CaloHitCollection *outputCollection = overlayer.GetOverlaidCollection();

		if( overlayer->OverlayDone() ) {

			const CaloHitCollection *outputCollection = overlayer->GetOverlaidCollection();
			// Fill the new LCIO CalorimeterHit collection
			IMPL::LCCollectionVec *lcOutputCollection = new IMPL::LCCollectionVec( LCIO::CALORIMETERHIT );
			UTIL::CellIDEncoder<IMPL::CalorimeterHitImpl> idEncoder (outputCodingPattern,lcOutputCollection);

			for( unsigned int i=0 ; i<outputCollection->size() ; i++ ) {

				CaloHit *caloHit = outputCollection->at(i);
				IMPL::CalorimeterHitImpl *hitImpl = new IMPL::CalorimeterHitImpl();

				float pos[3] = { float(caloHit->GetPosition().x())
								, float(caloHit->GetPosition().y())
								, float(caloHit->GetPosition().z()) };
				hitImpl->setPosition( pos );
				hitImpl->setType( caloHit->GetTypeID() );
				hitImpl->setTime( float(caloHit->GetTime()) );
				if( caloHit->GetThreshold() == fCaloHitThr1 ) hitImpl->setEnergy(1.0);
				if( caloHit->GetThreshold() == fCaloHitThr2 ) hitImpl->setEnergy(2.0);
				if( caloHit->GetThreshold() == fCaloHitThr3 ) hitImpl->setEnergy(3.0);
				idEncoder["I"] = caloHit->GetIJK().at(0);
				idEncoder["J"] = caloHit->GetIJK().at(1);
				idEncoder["K-1"] = caloHit->GetIJK().at(2);
				idEncoder.setCellID( hitImpl );

				lcOutputCollection->addElement( hitImpl );
			}

			// Fill the lcio track collection
			IMPL::LCCollectionVec *lcTrackOutputCollection = new IMPL::LCCollectionVec( LCIO::TRACK );

			if( overlayer->GetTrackPair().first != nullptr )
				lcTrackOutputCollection->addElement( overlayer->GetTrackPair().first );

			if( overlayer->GetTrackPair().second != nullptr )
				lcTrackOutputCollection->addElement( overlayer->GetTrackPair().second );

			// Fill the final event
			IMPL::LCEventImpl *outputEvent = new IMPL::LCEventImpl();

			IMPL::LCFlagImpl chFlag(0);
			EVENT::LCIO bitinfo;
			chFlag.setBit(bitinfo.CHBIT_LONG );   // calorimeter hit position
			chFlag.setBit(bitinfo.CHBIT_ID1 );    // cell ID
			chFlag.setBit(bitinfo.CHBIT_STEP );   // step info

			lcOutputCollection->setFlag( chFlag.getFlag()  ) ;
			nbOfOverlaidEvents++;

			outputEvent->addCollection( (LCCollection*) lcOutputCollection,outputCollectionName);
			outputEvent->addCollection( (LCCollection*) lcTrackOutputCollection , "SDHCALPrimaryTracks" );
			outputEvent->setRunNumber(0);
			outputEvent->setEventNumber(nbOfOverlaidEvents);
			outputEvent->setDetectorName( evt1->getDetectorName() );

			// write the final event
			lcWriter->writeEvent(outputEvent);

			sdhcal1->ClearContent();
			sdhcal2->ClearContent();

			nbOfLostHits += overlayer->GetNumberOfLostHits();
			evtID++;
		}
		else {

			sdhcal1->ClearContent();
			sdhcal2->ClearContent();
			nbOfSkippedEvents++;
			evtID++;
		}
	}


	/*******************
	 * Finalize the task
	 *******************/

	lcReader1->close();
	lcReader2->close();

	SdhcalConfig::Kill();

	delete parser;
	delete lcReader1;
	delete lcReader2;
	delete lcWriter;


	/*****************************
	 * Print the final statistics
	 *****************************/

	if(!discreteMode) {
		cout << endl;
		cout << "********************************************************"     << endl;
		cout << "*  OverlayEvent finished with the following statistics :"  << endl;
		cout << "*                                                       "  << endl;
		cout << "*  - Number of overlaid events : " << nbOfOverlaidEvents  << endl;
		cout << "*  - Number of skipped events   : " << nbOfSkippedEvents    << endl;
//		if(nbOfSkippedEvents != 0) {
//			cout << "*  - List of skipped events : ";
//			for(unsigned int i=0 ; i<skippedEventsVec.size() ; i++) {
//				cout << skippedEventsVec.at(i);
//				if(i!=skippedEventsVec.size()-1) cout << " , ";
//				else cout << endl;
//			}
//		}
		cout << "*  - Number of lost hits        : " << nbOfLostHits << endl;
		cout << "*  - Input sclio files          : " << endl;
		cout << "*       " << inputFile1 << endl;
		cout << "*       " << inputFile2 << endl;
		cout << "*  - Output sclio file          : " << endl;
		cout << "*       " << slcioOutputFile << endl;
		cout << "*  - Collection name            : " << outputCollectionName  << endl;
		cout << "*  - Separation distance        : " << separationDistance << " pads" << endl;
		cout << "*********************************************************" << endl;
	}

	return 0;
}

