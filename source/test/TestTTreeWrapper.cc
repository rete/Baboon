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
#include "Analysis/TTreeWrapper.hh"


using namespace std;
using namespace sdhcal;


int main (int argc ,char *argv[]) {

	cout << "Don't forget to source init_SDHCAL.sh script before running this..." << endl;


	/********************************************
	 * Grab the most useful environment variable
	 ********************************************/

//	string pathToSDHCAL = "";
//	pathToSDHCAL = getenv("PATH_TO_SDHCAL");
//	if( pathToSDHCAL.empty() ) throw runtime_error("'PATH_TO_SDHCAL' env variable is not set.\n Please source init_SDHCAL.sh before running.");


	TFile *f = new TFile("TestTTreeWrapper.root","RECREATE");

	TTreeWrapper *outputTTree = new TTreeWrapper();
//	treeWrapper->AddTree("HCALBarrel");
//	treeWrapper->AddTree("HCALBarrel_2");
//	treeWrapper->AddTree("HCALBarrel_3");



//	treeWrapper->AddBranch<int>("HCALBarrel","loopVar");
//	treeWrapper->AddBranch<int>("HCALBarrel_2","loopVar2");
//	treeWrapper->AddBranch<int>("HCALBarrel_3","loopVar3");


	cout << "debug" << endl;
	for(int i=0 ; i<200 ; i++) {

		cout << "debug i = " << i << endl;
		outputTTree->Set<int>("HCALBarrel","loopVar",i);
		outputTTree->Set<int>("HCALBarrel","loopVar2",i+1000);
		outputTTree->Set<int>("HCALBarrel_2","loopVar2",i+200);
		outputTTree->Set<int>("HCALBarrel_3","loopVar3",i+400);

		outputTTree->Fill("HCALBarrel");
		outputTTree->Fill("HCALBarrel_2");
		outputTTree->Fill("HCALBarrel_3");

	}

	cout << "debug before write" << endl;
	f->Write();
	cout << "debug before close" << endl;
	f->Close();
//	delete f;
	cout << "debug before delete wrapper" << endl;
//	delete treeWrapper;

	return EXIT_SUCCESS;
}
