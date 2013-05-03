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
#include "Config/SdhcalConfig.hh"

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
using namespace baboon;
//using namespace cfgparser;

// for debug messages (devel)
#define __DEBUG__ 0

enum HitThreshold {
	fThr1,
	fThr2,
	fThr3
};

HitThreshold GetHitThreshold(int thr);

int main (int argc ,char *argv[]) {

	/********************************************
	 * Grab the most useful environment variable
	 ********************************************/

	string pathToSDHCAL = getenv("PATH_TO_SDHCAL");
	if( pathToSDHCAL.empty() )
		throw runtime_error("'PATH_TO_SDHCAL' env variable is not set.\n "
		"Please source init_SDHCAL.sh before running.");

	string codingPattern ="";
	SdhcalConfig *config = SdhcalConfig::GetInstance();
	config->LoadFile( pathToSDHCAL + string("/config/SDHCAL.cfg") );
	config->GetData("general").GetValue("codingPattern",&codingPattern);


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
	TCLAP::ValueArg<std::string> collectionArg(
				"c"
				,"collection"
				,"collection name in slcio"
				,false
				,"HCALBarrel"
				,"string" );

	vector<string> allowedHitDisplay;
	allowedHitDisplay.push_back(string("threshold"));
	allowedHitDisplay.push_back(string("overlay"));
	TCLAP::ValuesConstraint<string> allowedValsHitDisplay( allowedHitDisplay );

	TCLAP::ValueArg<std::string> hitDisplayModeArg(
				""
				,"hit-display"
				,"Hit display mode"
				,false
				,"threshold"
				,&allowedValsHitDisplay);

	vector<string> allowedDisplay;
	allowedDisplay.push_back(string("3D"));
	allowedDisplay.push_back(string("XY"));
	TCLAP::ValuesConstraint<string> allowedValsDisplay( allowedDisplay );

	TCLAP::ValueArg<std::string> displayModeArg(
				""
				,"display-mode"
				,"Display mode for the third plot"
				,false
				,"XY"
				,&allowedValsDisplay);

	cmd.add( slcioFileArg );
	cmd.add( collectionArg );
	cmd.add( displayModeArg );
	cmd.add( hitDisplayModeArg );
	cmd.parse( argc, argv );

	TApplication app("Event Display of SDHCAL",&argc,argv);

	// Open the slcio file
	string slcioFile = slcioFileArg.getValue();
	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
	try{
		lcReader->open( slcioFile );
	} catch (IO::IOException& e ) {
		cerr << e.what() << endl;
		cerr << "DrawEvent have to stop. Sorry!" << endl;
		exit(1);
	}

	TCanvas *canvas = new TCanvas("canvas","Event display");
	canvas->SetCanvasSize(1100,650);
	canvas->SetWindowSize(1200,700);
	canvas->Divide(2,1);
	canvas->GetPad(1)->Divide(1,2);
	canvas->GetPad(2)->ResizePad();

	EVENT::LCEvent *event = NULL;
	string collectionName = collectionArg.getValue();


	while( ( event = lcReader->readNextEvent() )!=0 ) {

		EVENT::LCCollection *lcCollection = NULL;

		try {
			lcCollection = event->getCollection(collectionName);
		} catch (EVENT::DataNotAvailableException& e) {
			cerr << e.what() << endl;
			cerr << "DrawEvent have to stop. Sorry!" << endl;
			exit(1);
		}

		UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder(codingPattern);

		TH2 *xShowerProfileThr1 = new TH2D("xsp1","X shower profile",200,0,50,384,0,96);
		TH2 *xShowerProfileThr2 = new TH2D("xsp2","X shower profile",200,0,50,384,0,96);
		TH2 *xShowerProfileThr3 = new TH2D("xsp3","X shower profile",200,0,50,384,0,96);

		TH2 *yShowerProfileThr1 = new TH2D("ysp1","Y shower profile",200,0,50,384,0,96);
		TH2 *yShowerProfileThr2 = new TH2D("ysp2","Y shower profile",200,0,50,384,0,96);
		TH2 *yShowerProfileThr3 = new TH2D("ysp3","Y shower profile",200,0,50,384,0,96);

		TH2 *xyShowerProfileTh = new TH2D("xyspth","X-Y shower profile (weighted with threshold)",96,0,96,96,0,96);
//		TH2 *xyShowerProfile = new TH2D("xysp","X-Y shower profile (weighted without threshold)",96,0,96,96,0,96);

		TH3 *view3DThr1 = new TH3D("v3D1","3D view of SDHCAL event",200,0,50,384,0,96,384,0,96);
		TH3 *view3DThr2 = new TH3D("v3D2","3D view of SDHCAL event",200,0,50,384,0,96,384,0,96);
		TH3 *view3DThr3 = new TH3D("v3D3","3D view of SDHCAL event",200,0,50,384,0,96,384,0,96);

		xShowerProfileThr1->SetMarkerStyle(20);
		xShowerProfileThr1->SetMarkerSize(.4);
		xShowerProfileThr1->SetMarkerColor(3);

		xShowerProfileThr2->SetMarkerStyle(20);
		xShowerProfileThr2->SetMarkerSize(.4);
		xShowerProfileThr2->SetMarkerColor(4);

		xShowerProfileThr3->SetMarkerStyle(20);
		xShowerProfileThr3->SetMarkerSize(.4);
		xShowerProfileThr3->SetMarkerColor(2);

		yShowerProfileThr1->SetMarkerStyle(20);
		yShowerProfileThr1->SetMarkerSize(.4);
		yShowerProfileThr1->SetMarkerColor(3);

		yShowerProfileThr2->SetMarkerStyle(20);
		yShowerProfileThr2->SetMarkerSize(.4);
		yShowerProfileThr2->SetMarkerColor(4);

		yShowerProfileThr3->SetMarkerStyle(20);
		yShowerProfileThr3->SetMarkerSize(.4);
		yShowerProfileThr3->SetMarkerColor(2);

		xyShowerProfileTh->SetFillColor(2);
//		xyShowerProfile->SetFillColor(2);

		view3DThr1->SetMarkerStyle(20);
		view3DThr1->SetMarkerSize(.4);
		view3DThr1->SetMarkerColor(3);

		view3DThr2->SetMarkerStyle(20);
		view3DThr2->SetMarkerSize(.4);
		view3DThr2->SetMarkerColor(4);

		view3DThr3->SetMarkerStyle(20);
		view3DThr3->SetMarkerSize(.4);
		view3DThr3->SetMarkerColor(2);

		for( int i=0 ; i<lcCollection->getNumberOfElements() ; i++ ) {
			EVENT::CalorimeterHit *hit = static_cast<EVENT::CalorimeterHit*> ( lcCollection->getElementAt(i) );
			HitThreshold fThr = GetHitThreshold ( int(hit->getEnergy()) );
			int type = hit->getType();
			int I = idDecoder(hit)["I"];
			int J = idDecoder(hit)["J"];
			int K = idDecoder(hit)["K-1"];


			if( hitDisplayModeArg.getValue() == "overlay" ) {

				if(type == 1) {
					cout << "type 1" << endl;
					xShowerProfileThr1->Fill( K , I );
					yShowerProfileThr1->Fill( K , J );
					if( displayModeArg.getValue() == "3D" )
						view3DThr1->Fill( K , I , J );
					else if( displayModeArg.getValue() == "XY" )
						xyShowerProfileTh->Fill( I , J , 1 );
				}
				else if(type == 2) {
					cout << "type 2" << endl;
					xShowerProfileThr2->Fill( K , I );
					yShowerProfileThr2->Fill( K , J );
					if( displayModeArg.getValue() == "3D" )
						view3DThr2->Fill( K , I , J );
					else if( displayModeArg.getValue() == "XY" )
						xyShowerProfileTh->Fill( I , J , 5 );

				}
				else if(type == 3) {
					xShowerProfileThr3->Fill( K , I );
					yShowerProfileThr3->Fill( K , J );
					if( displayModeArg.getValue() == "3D" )
						view3DThr3->Fill( K , I , J );
					else if( displayModeArg.getValue() == "XY" )
						xyShowerProfileTh->Fill( I , J , 10 );
				}
			}


			if( hitDisplayModeArg.getValue() == "threshold" ) {

				if(fThr == fThr1) {
					xShowerProfileThr1->Fill( K , I );
					yShowerProfileThr1->Fill( K , J );
					if( displayModeArg.getValue() == "3D" )
						view3DThr1->Fill( K , I , J );
					else if( displayModeArg.getValue() == "XY" )
						xyShowerProfileTh->Fill( I , J , 1 );
				}
				else if(fThr == fThr2) {
					xShowerProfileThr2->Fill( K , I );
					yShowerProfileThr2->Fill( K , J );
					if( displayModeArg.getValue() == "3D" )
						view3DThr2->Fill( K , I , J );
					else if( displayModeArg.getValue() == "XY" )
						xyShowerProfileTh->Fill( I , J , 5 );
				}
				else if(fThr == fThr3) {
					xShowerProfileThr3->Fill( K , I );
					yShowerProfileThr3->Fill( K , J );
					if( displayModeArg.getValue() == "3D" )
						view3DThr3->Fill( K , I , J );
					else if( displayModeArg.getValue() == "XY" )
						xyShowerProfileTh->Fill( I , J , 10 );
				}
			}

		}

		gStyle->SetOptStat(0);

		canvas->GetPad(1)->cd(1);
		xShowerProfileThr1->Draw();
		xShowerProfileThr2->Draw("same");
		xShowerProfileThr3->Draw("same");

		canvas->GetPad(1)->cd(2);
		yShowerProfileThr1->Draw();
		yShowerProfileThr2->Draw("same");
		yShowerProfileThr3->Draw("same");

		canvas->cd(2);
		TView *viewer = new TView3D();
		viewer->SetRange(0,0,0,96,96,60);
		canvas->SetView(viewer);
		if( displayModeArg.getValue() == "XY" )
			xyShowerProfileTh->Draw("colz");
		else if( displayModeArg.getValue() == "3D" ) {
			view3DThr1->Draw();
			view3DThr2->Draw("same");
			view3DThr3->Draw("same");
		}

		ostringstream oss;
		oss << event->getEventNumber();
		string canvasTitle = "Event display : Event no " + oss.str() ;

		canvas->SetTitle( canvasTitle.c_str() );

		canvas->Update();

		canvas->WaitPrimitive();

		delete xShowerProfileThr1;
		delete xShowerProfileThr2;
		delete xShowerProfileThr3;

		delete yShowerProfileThr1;
		delete yShowerProfileThr2;
		delete yShowerProfileThr3;

//		delete xyShowerProfile;
		delete xyShowerProfileTh;

		delete view3DThr1;
		delete view3DThr2;
		delete view3DThr3;

	}

	SdhcalConfig::Kill();
//	app.Run();
	return 0;
}


HitThreshold GetHitThreshold(int thr) {
	if(thr == 2) return fThr1;
	if(thr == 1) return fThr2;
	if(thr == 3) return fThr3;
}





