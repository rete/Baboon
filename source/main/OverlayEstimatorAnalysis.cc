/*
 *
 * //OverlayEstimatorAnalysis.cc main file
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

// sdhcal includes
#include "Analysis/InputTreeWrapper.hh"
#include "Analysis/BasicDraw.hh"
#include "Analysis/OverlayEstimator/TreeProcessor.hh"

// root includes
#include "TCanvas.h"
#include "TGraph.hh"
#include "TTree.h"
#include "TFile.h"

// tclap includes
#include "tclap/CmdLine.h"

using namespace std;
using namespace cfgparser;
using namespace baboon;

// for debug messages (devel)
#define __DEBUG__ 0



int main (int argc ,char *argv[]) {

	cout << "Don't forget to source init_Baboon.sh script before running this..." << endl;
	/********************************************
	 * Grab the most useful environment variable
	 ********************************************/

	string pathToBaboon = "";
	pathToBaboon = getenv("BABOON_HOME");
	if( pathToBaboon.empty() ) throw runtime_error("'BABOON_HOME' env variable is not set.\n Please source init_SDHCAL.sh before running.");


	/*********************************
	 * Define the command line parser
	 *********************************/

	string footer = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine cmd(footer, ' ', "1.0");

	TCLAP::ValueArg<std::string> pathToRootFilesArg(
				"o"
				,"root-output"
				,"root output file"
				,true
				,"OverlayEstimator.root"
				,"string" );

	TCLAP::ValueArg<int> energyShower1Arg(
				"e"
				,"energy-shower"
				,"energy of one of the two showers"
				,true
				,""
				,"int" );

	cmd.add( pathToRootFilesArg );
	cmd.add( energyShower1Arg );
	cmd.parse( argc, argv );



	vector<int> energies;
	energies.push_back(10);
	energies.push_back(20);
	energies.push_back(30);
	energies.push_back(40);
	energies.push_back(50);

	vector<int> pads;
	pads.push_back(5);
	pads.push_back(10);
	pads.push_back(15);
	pads.push_back(20);
	pads.push_back(25);
	pads.push_back(30);

	TGraph* graphPadsPurity5 = NewTGraph(5,1);
	TGraph* graphPadsPurity10 = NewTGraph(5,2);
	TGraph* graphPadsPurity15 = NewTGraph(5,3);
	TGraph* graphPadsPurity20 = NewTGraph(5,4);
	TGraph* graphPadsPurity25 = NewTGraph(5,5);
	TGraph* graphPadsPurity30 = NewTGraph(5,6);

	TGraph* graphPadsContamination5 = NewTGraph(5,1);
	TGraph* graphPadsContamination10 = NewTGraph(5,2);
	TGraph* graphPadsContamination15 = NewTGraph(5,3);
	TGraph* graphPadsContamination20 = NewTGraph(5,4);
	TGraph* graphPadsContamination25 = NewTGraph(5,5);
	TGraph* graphPadsContamination30 = NewTGraph(5,6);



	for( unsigned int e=0 ; e<energies.size() ; e++ ) {

		for(  )


	}




	return 0;
}

