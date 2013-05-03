/*
 *
 * //Event selector main file
 * Copyright (c) 2012 CNRS / IPNL
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

// CfgParser includes
#include "Parser/CfgParser.hh"
#include "Parser/Data.hh"

// lcio includes
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCCollection.h"
#include "EVENT/SimCalorimeterHit.h"
#include "EVENT/LCObject.h"
#include "EVENT/Cluster.h"
#include "UTIL/LCTOOLS.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCTOOLS.h"

// tclap includes
#include "tclap/CmdLine.h"

// sdhcal includes
#include "Managers/AlgorithmManager.hh"
#include "Managers/AnalysisManager.hh"
#include "Algorithm/AbstractAlgorithm.hh"
#include "Algorithm/Tracking/HoughTransformAlgorithm.hh"
#include "Objects/Hit.hh"
#include "Managers/HitManager.hh"
#include "Objects/Cluster.hh"
#include "Managers/ClusteringManager.hh"
#include "Objects/TrackSegment.hh"
#include "SdhcalConfig.hh"

using namespace std;
using namespace sdhcal;
using namespace cfgparser;


int main (int argc ,char *argv[]) {

	cout << "Don't forget to source init_SDHCAL.sh script before running this..." << endl;


	/********************************************
	 * Grab the most useful environment variable
	 ********************************************/

	string pathToSDHCAL = "";
	pathToSDHCAL = getenv("PATH_TO_SDHCAL");
	if( pathToSDHCAL.empty() ) throw runtime_error("'PATH_TO_SDHCAL' env variable is not set.\n Please source init_SDHCAL.sh before running.");


	SdhcalConfig *config = SdhcalConfig::GetInstance();
	config->LoadFile( pathToSDHCAL + "/config/SDHCAL.cfg" );

	/*********************************
	 * Define the command line parser
	 *********************************/

	string footer = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine cmd(footer, ' ', "1.0");

	TCLAP::ValueArg<std::string> configFileArg(
				"f"
				,"config-file"
				,"config file to start the program"
				,true
				,""
				,"string" );

	cmd.add( configFileArg );
	cmd.parse( argc, argv );

	/****************************************
	 * Read the Event Selector configuration
	 *     of SDHCAL in EventSelector.cfg
	 ****************************************/

	string inputFileName = "";
	string outputFileName = "";
	string inputCollectionName = "";
	string outputCollectionName = "";

	CfgParser *parser = new CfgParser();
	parser->SetConfigFileName( configFileArg.getValue() );
	parser->Read();

	parser->GetValue( "input","slcioFile",&inputFileName );
	parser->GetValue( "input","collectionName",&inputCollectionName );
	parser->GetValue( "output","slcioFile",&outputFileName );
	parser->GetValue( "output","collectionName",&outputCollectionName );


	/***********************************************************
	 * Define all algorithms, add them in the algorithm manager
	 ***********************************************************/

	AlgorithmManager *algorithmManager = AlgorithmManager::GetInstance();
	algorithmManager->SetConfigFileName(pathToSDHCAL + "/config/Algorithm.cfg");

	// Add the Hough Transform Algorithm for track reconstruction within the sdhcal
	algorithmManager->RegisterAlgorithm( new HoughTransformAlgorithm() );

	// Initialize it
	algorithmManager->Initialize();





	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
	lcReader->open( inputFileName );

	IO::LCWriter* lcWriter = IOIMPL::LCFactory::getInstance()->createLCWriter();
	lcWriter->open( outputFileName , LCIO::WRITE_NEW );

	EVENT::LCEvent *evt;
	int evtID = 0;
	int nbOfKeptEvts = 0;
	int nbOfSkippedEvents = 0;


	while( (evt = lcReader->readNextEvent()) != 0 ) {

//		if( evt->getEventNumber() > 200 ) break;

		cout << "evt no " << evt->getEventNumber() << endl;

		bool keepEvent = false;

		EVENT::LCCollection *collection = evt->getCollection( inputCollectionName );

		HitManager *hitManager = HitManager::GetInstance();
		hitManager->BeginOfEvent(collection);
		HitCollection *hitCollection = hitManager->GetHitCollection();
		HitCollection *cutHitCollection = new HitCollection();

		for( unsigned int i=0 ; i<hitCollection->size() ;i++ ) {
			if( hitCollection->at(i)->GetIJK().at(2) <= 10 )
				cutHitCollection->push_back(hitCollection->at(i));
		}

		ClusteringManager *clusterManager = ClusteringManager::GetInstance();
		ClusterCollection *clusterCollection = clusterManager->GetCluster3D();

		AlgorithmManager *algoManager = AlgorithmManager::GetInstance();
		if( algoManager->AlgorithmIsRegistered( "HoughTransformAlgorithm" ) ) {

			HoughTransformAlgorithm *houghTransform = (HoughTransformAlgorithm *)algoManager->GetAlgorithm( "HoughTransformAlgorithm" );
			houghTransform->SetClusterCollection(clusterCollection);
			houghTransform->Process();
			TrackSegmentCollection *trackSegCollection = houghTransform->GetTrackSegmentCollection();

//			cout << "number of track found in the ten first layers : " << trackSegCollection->size() << endl;

			for( unsigned int i=0 ; i<trackSegCollection->size() ; i++ ) {

				TrackSegment *trackSegment = trackSegCollection->at(i);
				TrackSegmentExtremities extremties = trackSegment->GetExtremities();
				if( extremties.first->GetIJK().at(2) <= 2 || extremties.second->GetIJK().at(2) <= 2  ) {
					if( trackSegment->GetThetaAngle() > 2.8 )
						keepEvent = true;
				}
			}
		}


		if( !keepEvent ) {

			nbOfSkippedEvents++;
			delete cutHitCollection;
			evtID++;
//			cout << "  !!!! Event skipped !!!!  " << endl;
			continue;
		}

//		cout << "  !!!! Event kept !!!!  " << endl;
//		IMPL::LCCollectionVec *outputLCCollection = new IMPL::LCCollectionVec ( outputCollectionName );
//		IMPL::LCEventImpl *outputEvent = new IMPL::LCEventImpl();

//		for( unsigned int i=0 ; i<collection->getNumberOfElements() ;i++ ) {
//			IMPL::CalorimeterHitImpl *hit = static_cast<IMPL::CalorimeterHitImpl *> ( collection->getElementAt(i) );
//			outputLCCollection->addElement( hit );
////			cout << "hit x position : " << hit->getPosition()[0] << endl;
//		}

//		IMPL::LCFlagImpl chFlag(0) ;
//		EVENT::LCIO bitinfo;
//		chFlag.setBit(bitinfo.CHBIT_BARREL ) ;   // calorimeter hit position
//		chFlag.setBit(bitinfo.CHBIT_ID1 ) ;    // cell ID
//		chFlag.setBit(bitinfo.CHBIT_STEP ) ;   // step info
//
//		outputLCCollection->setFlag( chFlag.getFlag()  ) ;
//
//		outputEvent->addCollection( (EVENT::LCCollection*) outputLCCollection,outputCollectionName);
//		outputEvent->setRunNumber( evt->getRunNumber() );
//		outputEvent->setEventNumber(nbOfKeptEvts);
//		outputEvent->setDetectorName( evt->getDetectorName() );


		lcWriter->writeEvent( evt );


		evtID++;
		nbOfKeptEvts++;
		delete cutHitCollection;
		hitManager->EndOfEvent();
//		delete outputLCCollection;
//		delete outputEvent;

	}

	cout << endl;
	cout << "********************************************************"     << endl;
	cout << "*  EventSelector finished with the following statistics :"  << endl;
	cout << "*                                                       "  << endl;
	cout << "*  - Number of skipped events   : " << nbOfSkippedEvents    << endl;
	cout << "*  - Input sclio file           : " << endl;
	cout << "*       " << inputFileName << endl;
	cout << "*  - Output sclio file          : " << endl;
	cout << "*       " << outputFileName << endl;
	cout << "*  - Collection name            : " << inputCollectionName  << endl;
	cout << "*********************************************************" << endl;

	delete lcWriter;
	delete lcReader;
	SdhcalConfig::Kill();

	return EXIT_SUCCESS;
}
