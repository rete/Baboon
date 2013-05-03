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
#include "TGraph.h"
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

	TCLAP::ValueArg<std::string> cfgFileArg(
				"f"
				,"config"
				,"config file for analysis"
				,true
				,""
				,"string" );


	cmd.add( cfgFileArg );
	cmd.parse( argc, argv );


	CfgParser *parser = new CfgParser( cfgFileArg.getValue() );
	parser->Read();

	vector<string> energy1;
	vector<string> energy2;
	vector<string> pads;
	string pathToRootFiles;

	parser->GetValue("variables","energy1",&energy1);
	parser->GetValue("variables","energy2",&energy2);
	parser->GetValue("variables","pads",&pads);
	parser->GetValue("paths","pathToRootFiles",&pathToRootFiles);

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


	for( unsigned int e1=0 ; e1<energy1.size() ; e1++ ) {
		for( unsigned int e2=0 ; e2<energy2.size() ; e2++ ) {
			for( unsigned int p=0 ; p<pads.size() ; p++ ) {


				string rootFile = pathToRootFiles
								+ string("/")
								+ energy1.at(e1)
								+ string("GeV/")
								+ string("double_calorimeterhit_pi_")
								+ energy1.at(e1)
								+ string("_")
								+ energy2.at(e2)
								+ string("GeV_I0_1_overlay_")
								+ pads.at(p)
								+ string("pads_new.root");

				TFile *file = new TFile( rootFile.c_str() );
				TTree *tree = (TTree *) file->Get( rootFile.c_str() );

				TreeProcessor proc( tree );
				proc.Loop();
				EstimatorVars vars = proc.GetMeans();


			}
		}
	}




	return 0;
}

