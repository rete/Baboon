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
#include "UTIL/LCTOOLS.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCTOOLS.h"

// sdhcal includes
#include "Utils/Converters.hh"
#include "Utils/HitComputing.hh"
#include "Overlay/Overlayer.hh"



// root includes
#include "TCanvas.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

// tclap includes
#include "tclap/CmdLine.h"

using namespace std;
using namespace cfgparser;
using namespace baboon;

// for debug messages (devel)
#define __DEBUG__ 0


int main (int argc ,char *argv[]) {

//	cout << "Don't forget to source init_SDHCAL.sh script before running this..." << endl;
//	/********************************************
//	 * Grab the most useful environment variable
//	 ********************************************/
//
//	string pathToSDHCAL = "";
//	pathToSDHCAL = getenv("PATH_TO_SDHCAL");
//	if( pathToSDHCAL.empty() ) throw runtime_error("'PATH_TO_SDHCAL' env variable is not set.\n Please source init_SDHCAL.sh before running.");


	/*********************************
	 * Define the command line parser
	 *********************************/

	string footer = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine cmd(footer, ' ', "1.0");
	TCLAP::ValueArg<std::string> pathToSdhcalArg(
				"p"
				,"path-to-sdhcal"
				,"path to sdhcal to run the script"
				,false
				,"/home/remi/ilcsoft/SDHCAL"
				,"string" );
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
	cmd.add( pathToSdhcalArg );
	cmd.parse( argc, argv );

	string pathToSDHCAL = pathToSdhcalArg.getValue();
	bool discreteMode = discreteModeArg.getValue();


	/*********************************************************************
	 * Open a CfgParser and read the global configuration (SDHCAL.cfg)
	 * and the configuration file to overlay events. The user must
	 * provide the last one while the first one contains default settings
	 * about SDHCAL. Do not modify this one!
	 *********************************************************************/

	vector<int> nbOfPadsXYZ;
	SdhcalConfig *config = SdhcalConfig::GetInstance();
	config->LoadFile( pathToSDHCAL + "/config/SDHCAL.cfg" );
	config->GetData("pads").GetValue("nbOfPadsXYZ",&nbOfPadsXYZ);


	// get the input configuration file from command line argv[1]
	string inputFile1 = "";
	string inputCollectionName1 = "";
	string inputFile2 = "";
	string inputCollectionName2 = "";
	string slcioOutputFile = "";
	string outputCollectionName = "";
	int separationDistance = 0;
	int nbOfEventsToOverlay = 0;

	string overlayCfgFileName = "";

	overlayCfgFileName = cfgFileArg.getValue();


	CfgParser *parser = new CfgParser( overlayCfgFileName );
	parser->Read();

	parser->GetValue("input1","slcioFile",&inputFile1);
	parser->GetValue("input1","collectionName",&inputCollectionName1);

	parser->GetValue("input2","slcioFile",&inputFile2);
	parser->GetValue("input2","collectionName",&inputCollectionName2);

	parser->GetValue("output","nbOfEventsToOverlay",&nbOfEventsToOverlay);
	parser->GetValue("output","separationDistance",&separationDistance);
	parser->GetValue("output","slcioFile",&slcioOutputFile);
	parser->GetValue("output","collectionName",&outputCollectionName);


	/*****************************************
	 * LC readers (input) and writer (output)
	 *****************************************/

	IO::LCReader* lcReader1 = IOIMPL::LCFactory::getInstance()->createLCReader();
	IO::LCReader* lcReader2 = IOIMPL::LCFactory::getInstance()->createLCReader();
	IO::LCWriter* lcWriter  = IOIMPL::LCFactory::getInstance()->createLCWriter();
	lcReader1->open( inputFile1 );
	lcReader2->open( inputFile2 );
	lcWriter->open( slcioOutputFile , LCIO::WRITE_NEW );
	EVENT::LCEvent *evt1;
	EVENT::LCEvent *evt2;

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
	 * Loop over events. Overlap two events from two input files.
	 * If by redefining a hit position the hit is out of the sdhcal
	 * the hit is deleted form the output collection. The list of skipped
	 * is printed at the end for the user. Stop overlapping when the a
	 * given number of events are overlapped (see cfg file).
	 ********************************************************************/

	int nbOfSkippedEvents = 0;
	int nbOfOverlaidEvents = 0;
	int evtID=0;
	vector<int> skippedEventsVec;
	int lastPercentDisplayed = 0;
	int nbOfLostHits = 0;
	nbOfEventsToOverlay = min(maxNbOfEvents,nbOfEventsToOverlay);


	while( nbOfOverlaidEvents != nbOfEventsToOverlay ) {

		// avoid to be out of the collection range
//		if(evtID == maxNbOfEvents) break;

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


		/************************************************
		 * Make a copy of the input collections to avoid
		 * read-only exception while overlapping events.
		 * Each collection is tagged via hit->setType(1)
		 * or hit->setType(2)
		 ************************************************/

		EVENT::LCCollection *lcCollection1 = evt1->getCollection(inputCollectionName1);
		EVENT::LCCollection *lcCollection2 = evt2->getCollection(inputCollectionName2);

		IMPL::LCCollectionVec *newLCCollection1 = new IMPL::LCCollectionVec(inputCollectionName1);
		IMPL::LCCollectionVec *newLCCollection2 = new IMPL::LCCollectionVec(inputCollectionName2);

		// Copy collection 1
		for(int hitID=0 ; hitID<lcCollection1->getNumberOfElements() ; hitID++) {
			IMPL::CalorimeterHitImpl *hitImpl = static_cast<IMPL::CalorimeterHitImpl*> (lcCollection1->getElementAt(hitID) );
			IMPL::CalorimeterHitImpl *newHitImpl = Converter::CopyCalorimeterHitImpl(hitImpl);
			newHitImpl->setType(1);  // Tag it as hit collection 1
			newLCCollection1->addElement( newHitImpl );
		}
		// Copy collection 2
		for(int hitID=0 ; hitID<lcCollection2->getNumberOfElements() ; hitID++) {
			IMPL::CalorimeterHitImpl *hitImpl = static_cast<IMPL::CalorimeterHitImpl*> (lcCollection2->getElementAt(hitID) );
			IMPL::CalorimeterHitImpl *newHitImpl = Converter::CopyCalorimeterHitImpl(hitImpl);
			newHitImpl->setType(2);  // Tag it as hit collection 2
			newLCCollection2->addElement( newHitImpl );
		}

		// HitCollectionComputer is used to compute general things
		// for hit collections like cog, first hit layer...
//		HitCollectionComputer hitCollectionComputer (newLCCollection1);

//		ThreeVector cog1 = hitCollectionComputer.GetCenterOfGravity();
		ThreeVector cog1 = GetCenterOfGravity(newLCCollection1);
//		hitCollectionComputer.SetLCCollection(newLCCollection2);
//		ThreeVector cog2 = hitCollectionComputer.GetCenterOfGravity();
		ThreeVector cog2 = GetCenterOfGravity(newLCCollection2);

		Overlayer overlayer(newLCCollection1,newLCCollection2);


		// Separate the two events by the given separation distance
		// and center it at the center of the sdhcal in the x direction.
		overlayer.SetTranslations(
				new ThreeVector( nbOfPadsXYZ.at(0)/2.0 + separationDistance / 2.0 - cog1.x()
								 ,nbOfPadsXYZ.at(1)/2.0 - cog1.y()
								,0)
			   ,new ThreeVector( nbOfPadsXYZ.at(0)/2.0 - separationDistance / 2.0 - cog2.x()
					   	   	    ,nbOfPadsXYZ.at(1)/2.0 - cog2.y()
					   	   	    ,0) );

		overlayer.OverlayCollections();



		// Grab the final overlaid collection and fill more
		// information in the collection and in the event.
		IMPL::LCCollectionVec *outputLCCollection = overlayer.GetOverlaidCollection();

		int count = 0;
		for( unsigned int i=0 ; i<outputLCCollection->getNumberOfElements() ; i ++) {
			IMPL::CalorimeterHitImpl *hitImpl = static_cast<IMPL::CalorimeterHitImpl*> (outputLCCollection->getElementAt(i) );
			if( hitImpl->getType() != 3 ) {
				count ++;
//				cout << "hitImpl->getType() : " << hitImpl->getType() <<  endl;
			}
		}
//		cout << "count after overlay : " << count << endl;
//		cout << "overlayer.GetNumberOfLostHits() : " << overlayer.GetNumberOfLostHits() << endl;


		IMPL::LCEventImpl *outputEvent = new IMPL::LCEventImpl();

		IMPL::LCFlagImpl chFlag(0) ;
		EVENT::LCIO bitinfo;
		chFlag.setBit(bitinfo.CHBIT_LONG ) ;   // sim calorimeter hit position
		chFlag.setBit(bitinfo.CHBIT_ID1 ) ;    // cell ID
		chFlag.setBit(bitinfo.CHBIT_STEP ) ;   // step info

		outputLCCollection->setFlag( chFlag.getFlag()  ) ;
		nbOfOverlaidEvents++;

		outputEvent->addCollection( (LCCollection*) outputLCCollection,outputCollectionName);
		outputEvent->setRunNumber(0);
		outputEvent->setEventNumber(nbOfOverlaidEvents);
		outputEvent->setDetectorName( evt1->getDetectorName() );

		// write the final event
		lcWriter->writeEvent(outputEvent);

//		delete outputEvent;
		nbOfLostHits += overlayer.GetNumberOfLostHits();
		evtID++;
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
		if(nbOfSkippedEvents != 0) {
			cout << "*  - List of skipped events : ";
			for(unsigned int i=0 ; i<skippedEventsVec.size() ; i++) {
				cout << skippedEventsVec.at(i);
				if(i!=skippedEventsVec.size()-1) cout << " , ";
				else cout << endl;
			}
		}
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

