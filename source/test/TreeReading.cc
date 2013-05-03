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
#include <stdexcept>


// root includes
#include <TFile.h>
#include <TTree.h>

// sdhcal includes
#include "Managers/AnalysisManager.hh"
#include "Analysis/InputTreeWrapper.hh"


using namespace std;
using namespace sdhcal;


int main (int argc ,char *argv[]) {


	TFile *f = new TFile("TestTTreeWrapper.root");
	TTree *t = (TTree *) f->Get("HCALBarrel");

	InputTTreeWrapper *inputTTree = new InputTTreeWrapper(t);
	int nentries = inputTTree->GetNbOfEntries();
	int loop = 0;
	int loop2 = 0;
	for( int i=0 ; i<nentries ; i++ ) {

		inputTTree->LoadEntry(i);

		inputTTree->GetValue(string("loopVar"),loop);
		inputTTree->GetValue(string("loopVar2"),loop2);
		cout << "loop : " << loop << endl;
		cout << "loop2 : " << loop2 << endl;
	}

	f->Close();
	delete f;


	return EXIT_SUCCESS;
}
