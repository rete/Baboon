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

	vector<TObject*> storage;


	// Graphs and histograms declaration
	TH1D *measRecovSim5pads = NewTH1D( "measRecovSim5pads" , "" , 4*(energy1Val+energy2Val)+1 , -2*(energy1Val+energy2Val) , 2*(energy1Val+energy2Val) , 1 );
	TH1D *measRecovTB5pads  = NewTH1D( "measRecovTB5pads"  , "" , 4*(energy1Val+energy2Val)+1 , -2*(energy1Val+energy2Val) , 2*(energy1Val+energy2Val) , 2 );

	TH1D *measRecovSim30pads = NewTH1D( "measRecovSim30pads" , "" , 4*(energy1Val+energy2Val)+1 , -2*(energy1Val+energy2Val) , 2*(energy1Val+energy2Val) , 1 );
	TH1D *measRecovTB30pads  = NewTH1D( "measRecovTB30pads"  , "" , 4*(energy1Val+energy2Val)+1 , -2*(energy1Val+energy2Val) , 2*(energy1Val+energy2Val) , 2 );

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

	storage.push_back( measRecovSim5pads );
	storage.push_back( measRecovTB5pads );
	storage.push_back( measRecovSim30pads );
	storage.push_back( measRecovTB30pads );
	storage.push_back( meanRecMeasSim );
	storage.push_back( meanRecMeasTB );
	storage.push_back( sigmaRecMeasSim );
	storage.push_back( sigmaRecMeasTB );
	storage.push_back( sigma90RecMeasSim );
	storage.push_back( sigma90RecMeasTB );
	storage.push_back( proba2SigmaSim );
	storage.push_back( proba2SigmaTB );
	storage.push_back( proba3SigmaSim );
	storage.push_back( proba3SigmaTB );


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

			wrapperSim->GetValue( "showersFound" , showersFound );

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

	TCanvas *measRecov5padsCanvas = new TCanvas("measRecov5padsCanvas","Master Canvas for WaitPrimitive()");
	TCanvas *measRecov30padsCanvas = new TCanvas("measRecov30padsCanvas");
	TCanvas *meanRecMeasCanvas = new TCanvas("meanRecMeasCanvas");
	TCanvas *sigmaRecMeasCanvas = new TCanvas("sigmaRecMeasCanvas");
	TCanvas *sigma90RecMeasCanvas = new TCanvas("sigma90RecMeasCanvas");
	TCanvas *proba2SigmaCanvas = new TCanvas("proba2SigmaCanvas");
	TCanvas *proba3SigmaCanvas = new TCanvas("proba3SigmaCanvas");

	TLegend *measRecov5padsLegend = NewTLegend(0.55,0.65,0.76,0.82);
	TLegend *measRecov30padsLegend = NewTLegend(0.55,0.65,0.76,0.82);
	TLegend *meanRecMeasLegend = NewTLegend(0.55,0.65,0.76,0.82);
	TLegend *sigmaRecMeasLegend = NewTLegend(0.55,0.65,0.76,0.82);
	TLegend *sigma90RecMeasLegend = NewTLegend(0.55,0.65,0.76,0.82);
	TLegend *proba2SigmaLegend = NewTLegend(0.55,0.65,0.76,0.82);
	TLegend *proba3SigmaLegend = NewTLegend(0.55,0.65,0.76,0.82);

	TMultiGraph *meanRecMeasMulti = new TMultiGraph();
	TMultiGraph *sigmaRecMeasMulti = new TMultiGraph();
	TMultiGraph *sigma90RecMeasMulti = new TMultiGraph();
	TMultiGraph *proba2SigmaMulti = new TMultiGraph();
	TMultiGraph *proba3SigmaMulti = new TMultiGraph();

	DrawText("CALICE Preliminary");

	storage.push_back( measRecov5padsCanvas );
	storage.push_back( measRecov30padsCanvas );
	storage.push_back( meanRecMeasCanvas );
	storage.push_back( sigmaRecMeasCanvas );
	storage.push_back( sigma90RecMeasCanvas );
	storage.push_back( proba2SigmaCanvas );
	storage.push_back( proba3SigmaCanvas );

	storage.push_back( measRecov5padsLegend );
	storage.push_back( measRecov30padsLegend );
	storage.push_back( meanRecMeasLegend );
	storage.push_back( sigmaRecMeasLegend );
	storage.push_back( sigma90RecMeasLegend );
	storage.push_back( proba2SigmaLegend );
	storage.push_back( proba3SigmaLegend );

	storage.push_back( meanRecMeasMulti );
	storage.push_back( sigmaRecMeasMulti );
	storage.push_back( sigma90RecMeasMulti );
	storage.push_back( proba2SigmaMulti );
	storage.push_back( proba3SigmaMulti );


	string legendTitle = energy2 + "-GeV track";
	string testBeamLegend = "Data";

	CaliceStyle();

	measRecov5padsCanvas->cd();
	measRecov5padsLegend->AddEntry( measRecovSim5pads , G4PhysicsList.c_str() , "l" );
	measRecov5padsLegend->AddEntry( measRecovTB5pads , testBeamLegend.c_str() , "l" );
	measRecovSim5pads->Draw();
	measRecovTB5pads->Draw("same");
	storage.push_back( DrawText("CALICE Preliminary") );
	measRecov5padsLegend->Draw();
	measRecovSim5pads->SetXTitle("Energy [GeV]");
	measRecovSim5pads->SetTitle("Mean Energy (Measured-Recovered) , d=5cm");

	measRecov30padsCanvas->cd();
	measRecov30padsLegend->AddEntry( measRecovSim30pads , G4PhysicsList.c_str() , "l" );
	measRecov30padsLegend->AddEntry( measRecovTB30pads , testBeamLegend.c_str() , "l" );
	measRecovSim30pads->Draw();
	measRecovTB30pads->Draw("same");
	storage.push_back( DrawText("CALICE Preliminary") );
	measRecov30padsLegend->Draw();
	measRecovSim30pads->SetXTitle("Energy [GeV]");
	measRecovSim30pads->SetTitle("Mean Energy (Measured-Recovered) , d=30cm");

	meanRecMeasCanvas->cd();
	meanRecMeasMulti->Add( meanRecMeasSim );
	meanRecMeasMulti->Add( meanRecMeasTB );
	meanRecMeasLegend->AddEntry( meanRecMeasSim , G4PhysicsList.c_str() , "l" );
	meanRecMeasLegend->AddEntry( meanRecMeasTB , testBeamLegend.c_str() , "l" );
	meanRecMeasMulti->Draw("alp");
	meanRecMeasLegend->Draw();
	storage.push_back( DrawText("CALICE Preliminary") );
	meanRecMeasMulti->GetYaxis()->SetTitle("Mean Energy (Measured-Recovered)");
	meanRecMeasMulti->GetXaxis()->SetTitle("Distance between shower axis [cm]");

	sigmaRecMeasCanvas->cd();
	sigmaRecMeasMulti->Add( sigmaRecMeasSim );
	sigmaRecMeasMulti->Add( sigmaRecMeasTB );
	sigmaRecMeasLegend->AddEntry( sigmaRecMeasSim , G4PhysicsList.c_str() , "l" );
	sigmaRecMeasLegend->AddEntry( sigmaRecMeasTB , testBeamLegend.c_str() , "l" );
	sigmaRecMeasMulti->Draw("alp");
	sigmaRecMeasLegend->Draw();
	storage.push_back( DrawText("CALICE Preliminary") );
	sigmaRecMeasMulti->GetYaxis()->SetTitle("Deviation #sigma (Measured-Recovered)");
	sigmaRecMeasMulti->GetXaxis()->SetTitle("Distance between shower axis [cm]");

	sigma90RecMeasCanvas->cd();
	sigma90RecMeasMulti->Add( sigma90RecMeasSim );
	sigma90RecMeasMulti->Add( sigma90RecMeasTB );
	sigma90RecMeasLegend->AddEntry( sigma90RecMeasSim , G4PhysicsList.c_str() , "l" );
	sigma90RecMeasLegend->AddEntry( sigma90RecMeasTB , testBeamLegend.c_str() , "l" );
	sigma90RecMeasMulti->Draw("alp");
	sigma90RecMeasLegend->Draw();
	storage.push_back( DrawText("CALICE Preliminary") );
	sigma90RecMeasMulti->GetYaxis()->SetTitle("Deviation #sigma_{90} (Measured-Recovered)");
	sigma90RecMeasMulti->GetXaxis()->SetTitle("Distance between shower axis [cm]");

	proba2SigmaCanvas->cd();
	proba2SigmaMulti->Add( proba2SigmaSim );
	proba2SigmaMulti->Add( proba2SigmaTB );
	proba2SigmaLegend->AddEntry( proba2SigmaSim , G4PhysicsList.c_str() , "l" );
	proba2SigmaLegend->AddEntry( proba2SigmaTB , testBeamLegend.c_str() , "l" );
	proba2SigmaMulti->Draw("alp");
	proba2SigmaLegend->Draw();
	storage.push_back( DrawText("CALICE Preliminary") );
	proba2SigmaMulti->GetYaxis()->SetTitle("Probability of recovering within 2 #sigma");
	proba2SigmaMulti->GetXaxis()->SetTitle("Distance between shower axis [cm]");

	proba3SigmaCanvas->cd();
	proba3SigmaMulti->Add( proba3SigmaSim );
	proba3SigmaMulti->Add( proba3SigmaTB );
	proba3SigmaLegend->AddEntry( proba3SigmaSim , G4PhysicsList.c_str() , "l" );
	proba3SigmaLegend->AddEntry( proba3SigmaTB , testBeamLegend.c_str() , "l" );
	proba3SigmaMulti->Draw("alp");
	proba3SigmaLegend->Draw();
	storage.push_back( DrawText("CALICE Preliminary") );
	proba3SigmaMulti->GetYaxis()->SetTitle("Probability of recovering within 3 #sigma");
	proba3SigmaMulti->GetXaxis()->SetTitle("Distance between shower axis [cm]");

	measRecov5padsCanvas->Update();
	measRecov30padsCanvas->Update();
	meanRecMeasCanvas->Update();
	sigmaRecMeasCanvas->Update();
	sigma90RecMeasCanvas->Update();
	proba2SigmaCanvas->Update();
	proba3SigmaCanvas->Update();

	measRecov5padsCanvas->WaitPrimitive();

	for( unsigned int i=0 ; i<storage.size() ; i++ )
		delete storage.at(i);

	storage.clear();

	return 0;
}

