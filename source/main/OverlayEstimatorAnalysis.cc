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
#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TTree.h"
#include "TFile.h"
#include "TMultiGraph.h"
#include "TPad.h"
#include "TLegend.h"

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

	char *pathToBab = NULL;
	pathToBab = getenv("BABOON_HOME");
	if( pathToBab == NULL )
		throw runtime_error("'BABOON_HOME' env variable is not set.\n "
		"Please source init_Baboon.sh before running.");

	string pathToBaboon(pathToBab);

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


	TApplication app("appOfTheDeath",0,0);

	CfgParser *parser = new CfgParser( cfgFileArg.getValue() );
	parser->Read();

//	vector<string> energy1;
	string energy1;
	vector<string> energy2;
	vector<string> pads;
	string pathToRootFiles;

	parser->GetValue("variables","energy1",&energy1);
	parser->GetValue("variables","energy2",&energy2);
	parser->GetValue("variables","pads",&pads);
	parser->GetValue("paths","pathToRootFiles",&pathToRootFiles);

	map< string , TGraph* > purityGraphs;
	map< string , TGraph* > contaminationGraphs;
	map< string , TGraph* > algorithmEfficiencyGraphs;


	for( unsigned int i=0 ; i<pads.size() ; i++ ) {
		purityGraphs[ pads.at(i) ] = NewTGraph(5,i+1);
		contaminationGraphs[ pads.at(i) ] = NewTGraph(5,i+1);
		algorithmEfficiencyGraphs[ pads.at(i) ] = NewTGraph(5,i+1);
	}


	for( unsigned int e2=0 ; e2<energy2.size() ; e2++ ) {
		for( unsigned int p=0 ; p<pads.size() ; p++ ) {


			string rootFile = pathToRootFiles
							+ string("/")
							+ energy1
							+ string("GeV/")
							+ string("double_calorimeterhit_pi_")
							+ energy1
							+ string("_")
							+ energy2.at(e2)
							+ string("GeV_I0_1_overlay_")
							+ pads.at(p)
							+ string("pads_new.root");

			string treeName = string("double_calorimeterhit_pi_")
							+ energy1
							+ string("_")
							+ energy2.at(e2)
							+ string("GeV_I0_1_overlay_")
							+ pads.at(p)
							+ string("pads_new.root");

			TTree *tree = 0;
			TFile *file = new TFile( rootFile.c_str() );
			tree = (TTree *) file->Get( treeName.c_str() );
			cout << "root file : " << rootFile.c_str() << endl;
			cout << "tree name : " << treeName.c_str() << endl;
			TreeProcessor proc( tree );
			proc.Loop();
			EstimatorVars vars = proc.GetMeans();

			(purityGraphs[ pads.at(p) ])->SetPoint( e2 , atoi( energy2.at(e2).c_str() ) , vars.purity1 );
			(contaminationGraphs[ pads.at(p) ])->SetPoint( e2 , atoi( energy2.at(e2).c_str() ) , vars.contamination1 );
			(algorithmEfficiencyGraphs[ pads.at(p) ])->SetPoint( e2 , atoi( energy2.at(e2).c_str() ) , vars.algorithmEfficiency*100 );

			delete file;

		}
	}


	string ccPadsPurityTitle = string("Purity for a ")+ energy1 +string(" GeV track");
	string ccPadsContaminationTitle = string("Contamination for a ")+ energy1 +string(" GeV track");
	string ccPadsAlgorithmEfficiencyTitle = string("Algorithm efficiency for a ")+ energy1 +string(" GeV track");

	TCanvas *ccPadsPurity = new TCanvas( "ccPadsPurity" , ccPadsPurityTitle.c_str() );
	TCanvas *ccPadsContamination = new TCanvas( "ccPadsContamination" , ccPadsContaminationTitle.c_str() );
	TCanvas *ccPadsAlgorithmEfficiency = new TCanvas( "ccPadsAlgorithmEfficiency" , ccPadsAlgorithmEfficiencyTitle.c_str() );



	TMultiGraph *multiPadsPurity = new TMultiGraph();
	TLegend *legendPadsPurity = NewTLegend(0.55,0.65,0.76,0.82);

	TMultiGraph *multiPadsContamination = new TMultiGraph();
	TLegend *legendPadsContamination = NewTLegend(0.55,0.65,0.76,0.82);

	TMultiGraph *multiPadsAlgorithmEfficiency = new TMultiGraph();
	TLegend *legendPadsAlgorithmEfficiency = NewTLegend(0.55,0.65,0.76,0.82);


	for( unsigned int i=0 ; i<pads.size() ; i++ ) {

		string legend = pads.at(i) + string(" cm");

		legendPadsPurity->AddEntry( (purityGraphs[ pads.at(i) ]) , legend.c_str() , "L" );
		multiPadsPurity->Add( (purityGraphs[ pads.at(i) ]) );

		legendPadsContamination->AddEntry( (contaminationGraphs[ pads.at(i) ]) , legend.c_str() , "L" );
		multiPadsContamination->Add( (contaminationGraphs[ pads.at(i) ]) );

		legendPadsAlgorithmEfficiency->AddEntry( (algorithmEfficiencyGraphs[ pads.at(i) ]) , legend.c_str() , "L" );
		multiPadsAlgorithmEfficiency->Add( (algorithmEfficiencyGraphs[ pads.at(i) ]) );
	}


	// Draw purity graph
	ccPadsPurity->cd();
	ccPadsPurity->SetLogy(1);
	CaliceStyle();
	multiPadsPurity->Draw("alp");
	TPaveText* calicePreliminary = DrawText("CALICE Preliminary");
	legendPadsPurity->Draw();

	multiPadsPurity->GetXaxis()->SetTitle("E [GeV]");
	multiPadsPurity->GetYaxis()->SetTitle("Purity");
	multiPadsPurity->GetXaxis()->SetTitleSize(0.04);
	multiPadsPurity->GetYaxis()->SetTitleSize(0.04);
	multiPadsPurity->GetYaxis()->SetRangeUser(1,105);

	gPad->Modified();


	// Draw contamination graph
	ccPadsContamination->cd();
	ccPadsContamination->SetLogy(1);
	CaliceStyle();
	multiPadsContamination->Draw("alp");
	calicePreliminary->Draw();
	legendPadsContamination->Draw();

	multiPadsContamination->GetXaxis()->SetTitle("E [GeV]");
	multiPadsContamination->GetYaxis()->SetTitle("Contamination");
	multiPadsContamination->GetXaxis()->SetTitleSize(0.04);
	multiPadsContamination->GetYaxis()->SetTitleSize(0.04);
	multiPadsContamination->GetYaxis()->SetRangeUser(0.4,105);

	gPad->Modified();


	// Draw algorithm efficiency graph
	ccPadsAlgorithmEfficiency->cd();
//	ccPadsAlgorithmEfficiency->SetLogy(1);
	CaliceStyle();
	multiPadsAlgorithmEfficiency->Draw("alp");
	calicePreliminary->Draw();
	legendPadsAlgorithmEfficiency->Draw();

	multiPadsAlgorithmEfficiency->GetXaxis()->SetTitle("E [GeV]");
	multiPadsAlgorithmEfficiency->GetYaxis()->SetTitle("Algorithm Efficiency");
	multiPadsAlgorithmEfficiency->GetXaxis()->SetTitleSize(0.04);
	multiPadsAlgorithmEfficiency->GetYaxis()->SetTitleSize(0.04);
	multiPadsAlgorithmEfficiency->GetYaxis()->SetRangeUser(0.4,105);

	gPad->Modified();


	// Update graphs
	ccPadsPurity->Update();
	ccPadsContamination->Update();
	ccPadsAlgorithmEfficiency->Update();

	ccPadsPurity->WaitPrimitive();
	ccPadsContamination->WaitPrimitive();
	ccPadsAlgorithmEfficiency->WaitPrimitive();



	// Deletion area
	for( unsigned int i=0 ; i<pads.size() ; i++ ) {
		delete purityGraphs[ pads.at(i) ];
		delete contaminationGraphs[ pads.at(i) ];
	}
	delete legendPadsPurity;
	delete multiPadsPurity;
	delete ccPadsPurity;

	delete legendPadsContamination;
	delete multiPadsContamination;
	delete ccPadsContamination;

	delete legendPadsAlgorithmEfficiency;
	delete multiPadsAlgorithmEfficiency;
	delete ccPadsAlgorithmEfficiency;

	delete calicePreliminary;


	return 0;
}

