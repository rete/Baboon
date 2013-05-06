/*
 *
 * //LcioToRoot.cc main file
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

// lcio includes
#include "IO/LCReader.h"
#include "IO/LCWriter.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCCollection.h"
#include "EVENT/CalorimeterHit.h"
#include "UTIL/LCTOOLS.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCTOOLS.h"
#include "UTIL/CellIDDecoder.h"

// root includes
#include "TCanvas.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

using namespace std;

// for debug messages (devel)
#define __DEBUG__ 0

#ifdef __MAKECINT__
#pragma link C++ class vector<float> +;
#endif

int main( int argc , char *argv[] ) {


	if( argc < 2 || argc > 3 ) {
		cout << "LcioToRoot takes one (slcio file) or two (slcio file + root file) argument(s) " << endl;
		exit(1);
	}

	string lcioFileName = string(argv[1]);
	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
	lcReader->open( lcioFileName );
	EVENT::LCEvent *evt;

	vector<int> I;
	vector<int> J;
	vector<int> K;
	vector<double> x;
	vector<double> y;
	vector<double> z;
	vector<double> energy;
	vector<double> energyError;
	vector<int> type;
	vector<int> cellID0;
	vector<int> cellID1;
	vector<double> time;
	int event;

	string rootFileName;
	if( argc == 3 ) rootFileName = string(argv[2]);
	else rootFileName = lcioFileName.substr( 0 , lcioFileName.size() - 6  ) + ".root";
	TFile *rootFile = TFile::Open( rootFileName.c_str() , "RECREATE" );
	TTree *outputTree = new TTree("HCALBarrel","HCALBarrel tree from slcio file");

	outputTree->Branch("event",&event);
	outputTree->Branch("I","std::vector<int>",&I);
	outputTree->Branch("J","std::vector<int>",&J);
	outputTree->Branch("K","std::vector<int>",&K);
	outputTree->Branch("x","std::vector<double>",&x);
	outputTree->Branch("y","std::vector<double>",&y);
	outputTree->Branch("z","std::vector<double>",&z);
	outputTree->Branch("cellID0","std::vector<int>",&cellID0);
	outputTree->Branch("cellID1","std::vector<int>",&cellID1);
	outputTree->Branch("energy","std::vector<double>",&energy);
	outputTree->Branch("energyError","std::vector<double>",&energyError);
	outputTree->Branch("type","std::vector<int>",&type);
	outputTree->Branch("time","std::vector<double>",&time);

	UTIL::CellIDDecoder<EVENT::CalorimeterHit> idDecoder("M:3,S-1:3,I:9,J:9,K-1:6");



	while( (evt = lcReader->readNextEvent()) != 0 ) {

		EVENT::LCCollection *lcCol = evt->getCollection("HCALBarrel");

		event = evt->getEventNumber();
		cout << "event no " << event << endl;

		for( int i=0 ; i<lcCol->getNumberOfElements() ; i++ ) {

			EVENT::CalorimeterHit *hit = static_cast<EVENT::CalorimeterHit *> ( lcCol->getElementAt(i) );
			I.push_back( idDecoder(hit)["I"] );
			J.push_back( idDecoder(hit)["J"] );
			K.push_back( idDecoder(hit)["K-1"] );
			x.push_back( hit->getPosition()[0] );
			y.push_back( hit->getPosition()[1] );
			z.push_back( hit->getPosition()[2] );
			energy.push_back( hit->getEnergy() );
			energyError.push_back( hit->getEnergyError() );
			cellID0.push_back( hit->getCellID0() );
			cellID1.push_back( hit->getCellID1() );
			type.push_back( hit->getType() );
			time.push_back( hit->getTime() );

		}

		outputTree->Fill();

		I.clear();
		J.clear();
		K.clear();
		x.clear();
		y.clear();
		z.clear();
		energy.clear();
		energyError.clear();
		cellID0.clear();
		cellID1.clear();
		type.clear();
		time.clear();

	}

	rootFile->Write();
	rootFile->Close();

	delete rootFile;

	return 0;
}
