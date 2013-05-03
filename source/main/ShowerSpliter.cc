/*
 *
 * //Shower separation main file
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
#include "CfgParser/CfgParser.hh"
#include "CfgParser/Data.hh"

// lcio includes
#include "IO/LCReader.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCCollection.h"

// root includes
#include <TGraph2D.h>
#include <TFile.h>
#include <TTree.h>
#include <TApplication.h>

// tclap includes
#include "tclap/CmdLine.h"

// sdhcal includes
#include "Managers/AlgorithmManager.hh"
#include "Algorithm/AlgorithmHeaders.hh"
#include "Managers/JobManager.hh"
#include "Managers/AnalysisManager.hh"
#include "Objects/Hit.hh"
#include "Objects/Cluster.hh"
#include "Geometry/Cone.hh"
#include "Geometry/ThreeVector.hh"

using namespace std;
using namespace baboon;


int main (int argc ,char *argv[]) {

	cout << "Don't forget to source init_SDHCAL.sh script before running this..." << endl;


	/********************************************
	 * Grab the most useful environment variable
	 ********************************************/

	string pathToSDHCAL = "";
	pathToSDHCAL = getenv("PATH_TO_SDHCAL");
	if( pathToSDHCAL.empty() ) throw runtime_error("'PATH_TO_SDHCAL' env variable is not set.\n Please source init_SDHCAL.sh before running.");


	TApplication app("app",0,0);
	/*********************************
	 * Define the command line parser
	 *********************************/

	string footer = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine cmd(footer, ' ', "1.0");
	TCLAP::ValueArg<std::string> slcioFileArg(
				"f"
				,"slcio-file"
				,"slcio sdhcal event file to treat"
				,true
				,""
				,"string" );

	TCLAP::ValueArg<std::string> rootFileNameArg(
				"o"
				,"root-output"
				,"root output file name"
				,false
				,pathToSDHCAL + "/output/ShowerSpliter/ShowerSpliter.root"
				,"string" );

	TCLAP::SwitchArg discreteModeArg(
				"d"
				,"discrete"
				,"Set the discrete mode. No screen printing."
				,false);

	cmd.add( slcioFileArg );
	cmd.add( discreteModeArg );
	cmd.add( rootFileNameArg );
	cmd.parse( argc, argv );


	/********************************
	 * Read the general configuration
	 *   of SDHCAL in SDHCAL.cfg
	 ********************************/

	SdhcalConfig *config = SdhcalConfig::GetInstance();
	config->LoadFile( pathToSDHCAL + "/config/SDHCAL.cfg" );


	/***********************************************************
	 * Define all algorithms, add them in the algorithm manager
	 ***********************************************************/

	AlgorithmManager *algorithmManager = AlgorithmManager::GetInstance();
	algorithmManager->SetConfigFileName(pathToSDHCAL + "/config/Algorithm.cfg");

	// Add the Hough Transform Algorithm for track reconstruction within the sdhcal
	algorithmManager->RegisterAlgorithm( new HoughTransformAlgorithm() );

	// Add isolation tagging algorithm
	algorithmManager->RegisterAlgorithm( new IsolationTaggingAlgorithm() );

	// Add clustering (2D) algorithm
	algorithmManager->RegisterAlgorithm( new NearbyClustering2DAlgorithm() );

	// Add clustering (3D) algorithm
	algorithmManager->RegisterAlgorithm( new NearbyClustering3DAlgorithm() );

	// Add core finder algorithm
	algorithmManager->RegisterAlgorithm( new CoreFinderAlgorithm() );

	// Add cone beginning algorithm
	algorithmManager->RegisterAlgorithm( new ConeBeginningAlgorithm() );

	// Add pca algorithm
//	algorithmManager->RegisterAlgorithm( new PrincipalComponentAnalysis() );

	// Initialize it
	algorithmManager->Initialize();

	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
	lcReader->open( slcioFileArg.getValue() );

	EVENT::LCEvent *evt;

	JobManager *jobManager = JobManager::GetInstance();
	jobManager->Init();

	unsigned int nbOfEventsToProcess = 40;
//			lcReader->getNumberOfEvents();


	/****************************************************************
	 * Define the analysis manager and give it the configuration file
	 ****************************************************************/

	AnalysisManager *analysisManager = AnalysisManager::GetInstance();
	analysisManager->SetRootFileName( rootFileNameArg.getValue() );
	analysisManager->Init();

	while( (evt = lcReader->readNextEvent()) != 0 ) {

		cout << "evt no " << evt->getEventNumber() << endl;

		jobManager->ProcessEvent(evt);
		if(evt->getEventNumber() >= nbOfEventsToProcess ) break;

	}

	analysisManager->End();
	jobManager->End();

	delete lcReader;

	return EXIT_SUCCESS;
}
