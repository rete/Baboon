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

// sdhcal includes
#include "SdhcalConfig.hh"

// lcio includes
#include "IO/LCReader.h"
#include "Exceptions.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCCollection.h"
#include "EVENT/CalorimeterHit.h"
#include "UTIL/CellIDDecoder.h"

// root includes
#include "TApplication.h"
#include "TCanvas.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TStyle.h"
#include "TView3D.h"


// tclap includes
#include "tclap/CmdLine.h"
#include "tclap/ValuesConstraint.h"

using namespace std;

// for debug messages (devel)
#define __DEBUG__ 0


int main (int argc ,char *argv[]) {

	/********************************************
	 * Grab the most useful environment variable
	 ********************************************/

	string pathToSDHCAL = getenv("PATH_TO_SDHCAL");
	if( pathToSDHCAL.empty() )
		throw runtime_error("'PATH_TO_SDHCAL' env variable is not set.\n "
		"Please source init_SDHCAL.sh before running.");

//	SdhcalConfig *config = SdhcalConfig::GetInstance();
//	config->LoadFile( pathToSDHCAL + string("/config/SDHCAL.cfg") );


	/*********************************
	 * Define the command line parser
	 *********************************/

	string footer = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine cmd(footer, ' ', "1.0");
	TCLAP::ValueArg<std::string> slcioFileArg(
				"f"
				,"slcio-file"
				,"slcio sdhcal event file to display"
				,true
				,""
				,"string" );

	cmd.add( slcioFileArg );
	cmd.parse( argc, argv );


	// Open the slcio file
	string slcioFile = slcioFileArg.getValue();
	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
	try{
		lcReader->open( slcioFile );
	} catch (IO::IOException& e ) {
		cerr << e.what() << endl;
		cerr << "Executable have to stop. Sorry!" << endl;
		exit(1);
	}

	EVENT::LCEvent *event = 0;
	string collectionName = "HCALBarrel";

	const int nbOfEvents = lcReader->getNumberOfEvents();
	int evtID = 0;

	//	while( ( event = lcReader->readNextEvent() )!=0 ) {

	#pragma omp parallel for private(evtID,event) shared(lcReader,cout)
	for( evtID = 0 ; evtID<nbOfEvents ; evtID++ ) {

		#pragma omp critical
		{
			cout << "evtID : " << evtID << endl;
			event = lcReader->readNextEvent();
		}

		EVENT::LCCollection *lcCollection = event->getCollection("HCALBarrel");

		int size = lcCollection->getNumberOfElements();

		for( unsigned int eltID=0 ; eltID<size ; eltID++ ) {

			EVENT::CalorimeterHit *hit = static_cast<EVENT::CalorimeterHit*> ( lcCollection->getElementAt( eltID ) );
//			#pragma omp critical
//			{
//				cout << "prout" << endl;
//////				cout << "Hit at x=" << hit->getPosition()[0] << "  y=" << hit->getPosition()[1] << "  z=" << hit->getPosition()[2] << endl;
//////				cout << "  energy = " << hit->getEnergy() << endl;
//			}

		}






	}

//	SdhcalConfig::Kill();
	delete lcReader;
	return 0;
}




