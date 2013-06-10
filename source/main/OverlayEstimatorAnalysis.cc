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
#include "CfgParser/RawCfgParser.hh"
#include "CfgParser/Section.hh"

// baboon includes
#include "Analysis/InputTreeWrapper.hh"
#include "Analysis/BasicDraw.hh"
#include "Analysis/OverlayEstimator/TreeProcessor.hh"
#include "Utilities/Numeric.hh"

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

	// Mandatory for executables if graphs...
	TApplication app("appOfTheDeath",0,0);

	// Load the cfg file and grabs the variables
	string energy1;
	string energy2;
	vector<string> padDistances;
	string pathToTBRootFiles;
	string pathToSimRootFiles;
	string G4PhysicsList;
	int energy1Val = 0;
	int energy2Val = 0;
	vector<int> padDistancesVal;

	RawCfgParser *parser = new RawCfgParser();
	parser->Read( cfgFileArg.getValue() );
	parser->GetValue("paths","pathToTBRootFiles",&pathToTBRootFiles);
	parser->GetValue("paths","pathToSimRootFiles",&pathToSimRootFiles);
	parser->GetValue("variables","energy1",&energy1);
	parser->GetValue("variables","energy2",&energy2);
	parser->GetValue("variables","padDistances",&padDistances);
	parser->GetValue("variables","G4PhysicsList",&G4PhysicsList);

	energy1Val = atoi( energy1.c_str() );
	energy2Val = atoi( energy2.c_str() );
	for( unsigned int i=0 ; i<padDistances.size() ; i++ )
		padDistancesVal.push_back( atoi( padDistances.at(i).c_str() ) );


	// Graphs and histograms declaration
	TH1D *measRecovSim5pads = NewTH1D( "measRecovSim5pads" , "" , 2*max(energy1Val,energy2Val)+1 , -max(energy1Val,energy2Val) , max(energy1Val,energy2Val) , 1 );
	TH1D *measRecovTB5pads  = NewTH1D( "measRecovTB5pads"  , "" , 2*max(energy1Val,energy2Val)+1 , -max(energy1Val,energy2Val) , max(energy2Val,energy2Val) , 2 );

	TH1D *measRecovSim30pads = NewTH1D( "measRecovSim30pads" , "" , 2*max(energy1Val,energy2Val)+1 , -max(energy1Val,energy2Val) , max(energy1Val,energy2Val) , 1 );
	TH1D *measRecovTB30pads  = NewTH1D( "measRecovTB30pads"  , "" , 2*max(energy1Val,energy2Val)+1 , -max(energy1Val,energy2Val) , max(energy1Val,energy2Val) , 2 );

	TGraph *meanRecMeasSim = NewTGraph( padDistances.size() , 1 );
	TGraph *meanRecMeasTB  = NewTGraph( padDistances.size() , 2 );

	TGraph *sigmaRecMeasSim = NewTGraph( padDistances.size() , 1 );
	TGraph *sigmaRecMeasTB  = NewTGraph( padDistances.size() , 2 );

	TGraph *sigma90RecMeasSim = NewTGraph( padDistances.size() , 1 );
	TGraph *sigma90RecMeasTB  = NewTGraph( padDistances.size() , 2 );

	TGraph *proba2SigmaSim = NewTGraph( padDistances.size() , 1 );
	TGraph *proba2SigmaTB  = NewTGraph( padDistances.size() , 2 );

	TGraph *proba3SigmaSim = NewTGraph( padDistances.size() , 1 );
	TGraph *proba3SigmaTB  = NewTGraph( padDistances.size() , 2 );

	vector<double> deltaRecoEnergySim;
	vector<double> deltaRecoEnergyTB;


	for( unsigned int p=0 ; p<padDistances.size() ; p++ ) {

		ostringstream ss;
		ss << "/pi-_" << energy1 << "_" << energy2 << "GeV_" << padDistances.at( p ) << "pads.root";

		TFile *fileSim = TFile::Open( ( pathToSimRootFiles + ss.str() ).c_str() );
		TFile *fileTB  = TFile::Open( ( pathToTBRootFiles  + ss.str() ).c_str() );

		TTree *treeSim = (TTree *) fileSim->Get("EstimatorVariables");
		TTree *treeTB  = (TTree *)  fileTB->Get("EstimatorVariables");

		InputTTreeWrapper *wrapperSim = new InputTTreeWrapper( treeSim );
		InputTTreeWrapper *wrapperTB  = new InputTTreeWrapper( treeTB );

		int nbOfEntriesSim = wrapperSim->GetNbOfEntries();
		int nbOfEntriesTB  = wrapperTB->GetNbOfEntries();


		for( int entry=0 ; entry<nbOfEntriesSim ; entry++ ) {

			wrapperSim->LoadEntry( entry );

			double deltaRecoEnergy1 = 0;
			bool showersFound = false;

			wrapperTB->GetValue( "showersFound" , showersFound );

			if( !showersFound )
				continue;

			wrapperSim->GetValue( "deltaRecoEnergy1" , deltaRecoEnergy1 );

			if( padDistancesVal.at(p) == 5 )
				measRecovSim5pads->Fill( deltaRecoEnergy1 );
			else if( padDistancesVal.at(p) == 30 )
				measRecovSim30pads->Fill( deltaRecoEnergy1 );

			deltaRecoEnergySim.push_back( deltaRecoEnergy1 );
		}

		for( int entry=0 ; entry<nbOfEntriesTB ; entry++ ) {

			wrapperTB->LoadEntry( entry );

			double deltaRecoEnergy1 = 0;
			bool showersFound = false;

			wrapperTB->GetValue( "showersFound" , showersFound );

			if( !showersFound )
				continue;

			wrapperTB->GetValue( "deltaRecoEnergy1" , deltaRecoEnergy1 );

			if( padDistancesVal.at(p) == 5 )
				measRecovTB5pads->Fill( deltaRecoEnergy1 );
			else if( padDistancesVal.at(p) == 30 )
				measRecovTB30pads->Fill( deltaRecoEnergy1 );

			deltaRecoEnergyTB.push_back( deltaRecoEnergy1 );
		}

		meanRecMeasSim->SetPoint( p , padDistancesVal.at(p) , Mean( deltaRecoEnergySim ) );
		meanRecMeasTB->SetPoint( p , padDistancesVal.at(p) , Mean( deltaRecoEnergyTB ) );

		sigmaRecMeasSim->SetPoint( p , padDistancesVal.at(p) , RMS( deltaRecoEnergySim ) );
		sigmaRecMeasTB->SetPoint( p , padDistancesVal.at(p) , RMS( deltaRecoEnergyTB ) );

		sigma90RecMeasSim->SetPoint( p , padDistancesVal.at(p) , RMS90( deltaRecoEnergySim ) );
		sigma90RecMeasTB->SetPoint( p , padDistancesVal.at(p) , RMS90( deltaRecoEnergyTB ) );

		proba2SigmaSim->SetPoint( p , padDistancesVal.at(p) , RecoveryProbabilityWithinSigma( deltaRecoEnergySim , 2 ) );
		proba2SigmaTB->SetPoint( p , padDistancesVal.at(p) , RecoveryProbabilityWithinSigma( deltaRecoEnergyTB , 2 ) );

		proba3SigmaSim->SetPoint( p , padDistancesVal.at(p) , RecoveryProbabilityWithinSigma( deltaRecoEnergySim , 3 ) );
		proba3SigmaTB->SetPoint( p , padDistancesVal.at(p) , RecoveryProbabilityWithinSigma( deltaRecoEnergyTB , 3 ) );

		fileSim->Close();
		fileTB->Close();
		delete fileSim;
		delete fileTB;
	}

	vector<TLegend*> legendStorage;

	TCanvas *measRecov5padsCanvas = new TCanvas("measRecov5padsCanvas","Master Canvas for WaitPrimitive()");
	TCanvas *measRecov30padsCanvas = new TCanvas("measRecov30padsCanvas");
	TCanvas *meanRecMeasCanvas = new TCanvas("meanRecMeasCanvas");
	TCanvas *sigmaRecMeasCanvas = new TCanvas("sigmaRecMeasCanvas");
	TCanvas *sigma90RecMeasCanvas = new TCanvas("sigma90RecMeasCanvas");
	TCanvas *proba2SigmaCanvas = new TCanvas("proba2SigmaCanvas");
	TCanvas *proba3SigmaCanvas = new TCanvas("proba3SigmaCanvas");

	string legendTitle = energy2 + "-GeV track";
	string testBeamLegend = "Data";

	measRecov5padsCanvas->cd();
	measRecovSim5pads->Draw();
	measRecovTB5pads->Draw("same");

	measRecov30padsCanvas->cd();
	measRecovSim30pads->Draw();
	measRecovTB30pads->Draw("same");

	meanRecMeasCanvas->cd();
	meanRecMeasSim->Draw("p");
	meanRecMeasTB->Draw("same p");

	sigmaRecMeasCanvas->cd();
	sigmaRecMeasSim->Draw("p");
	sigmaRecMeasTB->Draw("same p");

	sigma90RecMeasCanvas->cd();
	sigma90RecMeasSim->Draw("p");
	sigma90RecMeasTB->Draw("same p");

	proba2SigmaCanvas->cd();
	proba2SigmaSim->Draw("p");
	proba2SigmaTB->Draw("same p");

	proba3SigmaCanvas->cd();
	proba3SigmaSim->Draw("p");
	proba3SigmaTB->Draw("same p");
	TLegend *proba3SigmaLegend = NewTLegend(0.55,0.65,0.76,0.82);
	proba3SigmaLegend->AddEntry( proba3SigmaSim , G4PhysicsList.c_str() , "L" );
	proba3SigmaLegend->AddEntry( proba3SigmaTB , testBeamLegend.c_str() , "L" );
	legendStorage.push_back( proba3SigmaLegend );

	measRecov5padsCanvas->Update();
	measRecov30padsCanvas->Update();
	meanRecMeasCanvas->Update();
	sigmaRecMeasCanvas->Update();
	sigma90RecMeasCanvas->Update();
	proba2SigmaCanvas->Update();
	proba3SigmaCanvas->Update();

	measRecov5padsCanvas->WaitPrimitive();

	for( unsigned int i=0 ; i<legendStorage.size() ; i++ )
		delete legendStorage.at(i);
	legendStorage.clear();

	delete measRecov5padsCanvas;
	delete measRecov30padsCanvas;
	delete meanRecMeasCanvas;
	delete sigmaRecMeasCanvas;
	delete sigma90RecMeasCanvas;
	delete proba2SigmaCanvas;
	delete proba3SigmaCanvas;

	/*
//	vector<string> energy1;
	string energy1;
	string energy2;
	vector<string> pads;


	parser->GetValue("variables","energy1",&energy1);
	parser->GetValue("variables","energy2",&energy2);
	parser->GetValue("variables","pads",&pads);
	parser->GetValue("paths","pathToTBRootFiles",&pathToTBRootFiles);
	parser->GetValue("paths","pathToSimRootFiles",&pathToSimRootFiles);

	map< string , TGraph* > shower1DeltaRecoveryGraphs;
	map< string , TGraph* > shower2DeltaRecoveryGraphs;
	map< string , TGraph* > algorithmEfficiencyGraphs;


	for( unsigned int i=0 ; i<pads.size() ; i++ ) {

		shower1DeltaRecoveryGraphs[ pads.at(i) ] = NewTGraph( 6 , i+1 );
		shower2DeltaRecoveryGraphs[ pads.at(i) ] = NewTGraph( 6 , i+1 );
		algorithmEfficiencyGraphs[ pads.at(i) ] = NewTGraph( 6 , i+1 );
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
*/

	return 0;
}

