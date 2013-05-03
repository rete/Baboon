/*
 *
 * //OverlayEstimator.cc main file
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

// lcio includes
#include "IO/LCReader.h"
#include "IO/LCWriter.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCCollection.h"
#include "UTIL/LCTOOLS.h"
#include "IMPL/LCCollectionVec.h"
#include "IMPL/LCEventImpl.h"
#include "IMPL/LCFlagImpl.h"
#include "IMPL/LCTOOLS.h"

// sdhcal includes
#include "Objects/Cluster.hh"
#include "Objects/Hit.hh"
#include "Utils/HitComputing.hh"
#include "Managers/HitManager.hh"
#include "Managers/AnalysisManager.hh"
#include "Geometry/Cylinder.hh"
#include "Managers/AlgorithmManager.hh"
#include "Algorithm/PrincipalComponentAnalysis.hh"
#include "Algorithm/IsolationTaggingAlgorithm.hh"


// root includes
#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TBranch.h"

// tclap includes
#include "tclap/CmdLine.h"

using namespace std;
using namespace cfgparser;
using namespace sdhcal;

// for debug messages (devel)
#define __DEBUG__ 0


//extern int gErrorIgnoreLevel;
//gErrorIgnoreLevel = kWarning;

int main (int argc ,char *argv[]) {

	cout << "Don't forget to source init_SDHCAL.sh script before running this..." << endl;
	/********************************************
	 * Grab the most useful environment variable
	 ********************************************/

	string pathToSDHCAL = "";
	pathToSDHCAL = getenv("PATH_TO_SDHCAL");
	if( pathToSDHCAL.empty() ) throw runtime_error("'PATH_TO_SDHCAL' env variable is not set.\n Please source init_SDHCAL.sh before running.");


	/*********************************
	 * Define the command line parser
	 *********************************/

	string footer = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine cmd(footer, ' ', "1.0");

	TCLAP::ValueArg<std::string> slcioFileArg(
				"f"
				,"slcio"
				,"slcio to run OverlayEstimator"
				,true
				,""
				,"string" );
	TCLAP::ValueArg<std::string> rootOutputFileArg(
				"o"
				,"root-output"
				,"root output file"
				,false
				,"OverlayEstimator.root"
				,"string" );

	TCLAP::SwitchArg gridModeArg(
				"g"
				,"grid"
				,"Set the grid mode. Allow to run OverlayEvent on the grid without using cin checking ."
				,false);

	vector<string> estimatorMode;
	estimatorMode.push_back(string("cylinder"));
	estimatorMode.push_back(string("pca"));
	TCLAP::ValuesConstraint<string> allowedestimatorMode( estimatorMode );

	TCLAP::ValueArg<std::string> estimatorModeArg(
				"m"
				,"estimator-mode"
				,"Estimator mode"
				,false
				,"cylinder"
				,&allowedestimatorMode);

	cmd.add( slcioFileArg );
	cmd.add( estimatorModeArg );
	cmd.add( gridModeArg );
	cmd.add( rootOutputFileArg );
	cmd.parse( argc, argv );

//	gErrorIgnoreLevel = kWarning;
//	TApplication app("app",0,0);


	vector<int> nbOfPadsXYZ;
	SdhcalConfig *config = SdhcalConfig::GetInstance();
	config->LoadFile( pathToSDHCAL + "/config/SDHCAL.cfg" );
	config->GetData("pads").GetValue("nbOfPadsXYZ",&nbOfPadsXYZ);


	/*************
	 * LC reader
	 *************/

	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
	lcReader->open( slcioFileArg.getValue() );
	EVENT::LCEvent *evt;

	/**********************************
	 * Algorithm manager instance and
	 * algorithm registrations.
	 **********************************/

	AlgorithmManager *algorithmManager = AlgorithmManager::GetInstance();
	algorithmManager->SetConfigFileName(pathToSDHCAL + "/config/Algorithm.cfg");

	// Add pca algorithm
	algorithmManager->RegisterAlgorithm( new PrincipalComponentAnalysis() );

	// Add isolation algorithm
	algorithmManager->RegisterAlgorithm( new IsolationTaggingAlgorithm() );

	// Initialize it
	algorithmManager->Initialize();


	/****************************
	 * Analysis settings
	 ****************************/

	AnalysisManager *analysisManager = AnalysisManager::GetInstance();
	analysisManager->SetRootFileName( rootOutputFileArg.getValue() );
	analysisManager->Init();

	/********************
	 * Loop over events.
	 ********************/

	int evtID = 0;
	int lastPercentDisplayed = 0;
	int nbOfEvts = lcReader->getNumberOfEvents();

//	TCanvas *cc = new TCanvas("c1","PCA");
//	cc->Divide(3,0);


	while( (evt = lcReader->readNextEvent()) != 0 ) {

		evtID = evt->getEventNumber();

		double overlaidHits = 0;

		// Progress display
		double percent = (evtID*1.0) / (1.0*nbOfEvts) * 10;
		if( round(percent) > lastPercentDisplayed ) {
			lastPercentDisplayed = round(percent);
			cout << "Processing ... " << lastPercentDisplayed*10 << "%" << endl;
		}


		HitManager *hitManager = HitManager::GetInstance();
		hitManager->BeginOfEvent( evt->getCollection("HCALBarrel") );

		HitCollection *hitCollection = hitManager->GetHitCollection();

		HitCollection *hitColType1 = new HitCollection();
		HitCollection *hitColType2 = new HitCollection();

		for( unsigned int i=0 ; i<hitCollection->size() ; i++ ) {

			if( hitCollection->at(i)->getType() == 1 ) {
//				cout << "type 1" << endl;
				hitColType1->push_back(hitCollection->at(i));
			}
			else if( hitCollection->at(i)->getType() == 2 ) {
//				cout << "type 2" << endl;
				hitColType2->push_back(hitCollection->at(i));
			}
			else if( hitCollection->at(i)->getType() == 3) {
				hitColType1->push_back(hitCollection->at(i));
				hitColType2->push_back(hitCollection->at(i));
				overlaidHits++;
			}
			else {
				cerr << "No shower type. Please check if the input slcio "
					"file has been produced by OverlayEvent script..." << endl;
				exit(1);
			}
		}


		double contamination1 = 0;
		double contamination2 = 0;

		double purity1 = 0;
		double purity2 = 0;

		double compensation1 = 0;
		double compensation2 = 0;

		double hitsFrom1in1 = 0;
		double hitsFrom1in2 = 0;
		double hitsFrom2in1 = 0;
		double hitsFrom2in2 = 0;

		bool showersFound = true;


		if( estimatorModeArg.getValue() == "pca" ) {

			/***************
			 * PC Analysis
			 ***************/

			TVectorD pcaEigenValues(2);
			TMatrixD pcaEigenVectors(2,2);

			if( algorithmManager->AlgorithmIsRegistered("PrincipalComponentAnalysis") ) {

				PrincipalComponentAnalysis *pca = (PrincipalComponentAnalysis *) algorithmManager->GetAlgorithm("PrincipalComponentAnalysis");
				pca->SetHitCollection( hitCollection );
				pca->Process();
				pcaEigenValues = pca->GetEigenValues();
				pcaEigenVectors = pca->GetEigenVectors();
			}
			else continue;

			ThreeVector cog = GetCenterOfGravity(hitCollection);
			TMatrixD initialDataSet(2,hitCollection->size());

			double *x = new double[ hitCollection->size() ];
			double *y = new double[ hitCollection->size() ];
			double *xNew = new double[ hitCollection->size() ];
			double *yNew = new double[ hitCollection->size() ];

			for( unsigned int i=0 ; i<hitCollection->size() ; i++  ) {
				initialDataSet(0,i) = hitCollection->at(i)->GetPosition().x() - cog.x();
				initialDataSet(1,i) = hitCollection->at(i)->GetPosition().y() - cog.y();
				x[i] = initialDataSet(0,i);
				y[i] = initialDataSet(1,i);
			}

			TMatrixD newDataSet(pcaEigenVectors,TMatrixD::kTransposeMult,initialDataSet);

			double xNewMean = 0;
			double yNewMean = 0;

			for( int i=0 ; i<hitCollection->size() ; i++) {
				xNewMean += newDataSet(0,i) / hitCollection->size();
				yNewMean += newDataSet(1,i) / hitCollection->size();
			}

			TH1D* newXHist = new TH1D("newXHist","newXHist",350,-700,700);

			for( unsigned int i=0 ; i<hitCollection->size() ; i++ ) {

				newDataSet(0,i) -= xNewMean;
				newDataSet(1,i) -= yNewMean;
				xNew[i] = newDataSet(0,i);
				yNew[i] = newDataSet(1,i);
				newXHist->Fill(newDataSet(0,i));
			}


			TGraph *dataBeforePCA = new TGraph(hitCollection->size(),x,y);
			dataBeforePCA->SetMarkerStyle(21);
			dataBeforePCA->SetMarkerSize(0.5);
			dataBeforePCA->SetMarkerColor(1);
			dataBeforePCA->SetLineWidth(2);
			dataBeforePCA->SetLineColor(1);

			TGraph *dataAfterPCA = new TGraph(hitCollection->size(),xNew,yNew);
			dataAfterPCA->SetMarkerStyle(21);
			dataAfterPCA->SetMarkerSize(0.5);
			dataAfterPCA->SetMarkerColor(2);
			dataAfterPCA->SetLineWidth(2);
			dataAfterPCA->SetLineColor(2);

//			cc->cd(3);
			newXHist->Smooth(400);
			TSpectrum *s1 = new TSpectrum( 2 );
			int nbOfPeaks = s1->Search(newXHist,2,"nodraw");

			if( nbOfPeaks <= 1 ) {
				showersFound = false;
				delete x;
				delete y;
				delete xNew;
				delete yNew;
				delete newXHist;
				delete s1;
//				continue;
			}
			else {


				double firstMean = min( s1->GetPositionX()[0] , s1->GetPositionX()[1]);
				double secondMean = max( s1->GetPositionX()[0] , s1->GetPositionX()[1]);
				double center = (secondMean - firstMean ) / 2.0 + firstMean ;

				double firstPosition = 0;
				double secondPosition = 0;
				double firstRMS = 0;
				double secondRMS = 0;
				int nbOfPosFirst = 0;
				int nbOfPosSecond = 0;

				double infLimitFirst = firstMean - abs( center - firstMean );
				double supLimitFirst = center;
				double infLimitSecond = center;
				double supLimitSecond = secondMean + abs( center - secondMean );

				for( unsigned int i=0 ; i<hitCollection->size() ; i++ ) {

					if(    (newDataSet(0,i) < supLimitFirst)
						&&  (newDataSet(0,i) > infLimitFirst  )       ) {

						firstRMS += (newDataSet(0,i) - firstMean)*(newDataSet(0,i) - firstMean);
						nbOfPosFirst ++;

					}

					else if(    (newDataSet(0,i) > infLimitSecond )
						&&  (newDataSet(0,i) < supLimitSecond)     ) {

						secondRMS += (newDataSet(0,i) - secondMean)*(newDataSet(0,i) - secondMean);
						nbOfPosSecond++;

					}
				}

				firstRMS = sqrt(firstRMS/ nbOfPosFirst);
				secondRMS = sqrt(secondRMS/ nbOfPosSecond);

				double realCenter = ( firstMean/firstRMS + secondMean/secondRMS ) / ( 1/firstRMS + 1/secondRMS );

				double commonHits = 0;

				double *xNew1 = new double[ hitColType1->size() ];
				double *xNew2 = new double[ hitColType2->size() ];
				double *yNew1 = new double[ hitColType1->size() ];
				double *yNew2 = new double[ hitColType2->size() ];
				int id1 = 0;
				int id2 = 0;

				for( unsigned int i=0 ; i<hitCollection->size() ; i++ ) {

					int type = hitCollection->at(i)->getType();

		//			if( type == 1 || type == 3 ) {
		//				xNew1[id1] = newDataSet(0,i);
		//				yNew1[id1] = newDataSet(1,i);
		//				id1++;
		//			}
		//			if( type == 2 || type == 3 ) {
		//				xNew2[id2] = newDataSet(0,i);
		//				yNew2[id2] = newDataSet(1,i);
		//				id2++;
		//			}

					if( newDataSet(0,i) < realCenter && (type == 1 || type == 3)) hitsFrom1in1++;
					if( newDataSet(0,i) < realCenter && (type == 2 || type == 3)) hitsFrom2in1++;
					if( newDataSet(0,i) > realCenter && (type == 2 || type == 3)) hitsFrom2in2++;
					if( newDataSet(0,i) > realCenter && (type == 1 || type == 3)) hitsFrom1in2++;
				}

				if( hitsFrom1in2 > hitsFrom2in2 ) {

					contamination1 = hitsFrom2in2 / (hitsFrom1in2 + hitsFrom2in2);
					contamination2 = hitsFrom1in1 / (hitsFrom2in1 + hitsFrom1in1);

					purity1 = ( hitsFrom1in2 ) / ( hitsFrom2in2 + hitsFrom1in2 );
					purity2 = ( hitsFrom2in1 ) / ( hitsFrom1in1 + hitsFrom2in1 );

					compensation1 = ( hitsFrom2in1 ) / hitColType2->size();
					compensation2 = ( hitsFrom1in2 ) / hitColType1->size();
				}
				else {
					contamination1 = hitsFrom2in1 / (hitsFrom2in1 + hitsFrom1in1);
					contamination2 = hitsFrom1in2 / (hitsFrom2in1 + hitsFrom2in2);

					purity1 = ( hitsFrom1in1 ) / ( hitsFrom2in1 + hitsFrom1in1 );
					purity2 = ( hitsFrom2in2 ) / ( hitsFrom1in2 + hitsFrom2in2 );

					compensation2 = ( hitsFrom2in1 ) / hitColType2->size();
					compensation1 = ( hitsFrom1in2 ) / hitColType1->size();
				}

		//		TMultiGraph multi;
		//		TGraph graph1( hitColType1->size(),xNew1,yNew1 );
		//		TGraph graph2( hitColType2->size(),xNew2,yNew2 );
		//
		//		graph1.SetMarkerStyle(21);
		//		graph1.SetMarkerSize(0.5);
		//		graph1.SetMarkerColor(2);
		//
		//		graph2.SetMarkerStyle(21);
		//		graph2.SetMarkerSize(0.5);
		//		graph2.SetMarkerColor(3);
		//
		//		multi.Add(&graph1);
		//		multi.Add(&graph2);


		//		TLine separationLine( realCenter , 0 , realCenter , 20 );
		//		TLine sep2(realCenter , -200 , realCenter , 200);
		//		separationLine.Draw();

		//		cc->cd(1);
		//		dataBeforePCA->Draw("ap");
		//
		//		cc->cd(2);
		////		dataAfterPCA->Draw("ap");
		//		multi.Draw("ap");
		//		sep2.Draw();

		//		cout << "hitsFrom1in1 : " << hitsFrom1in1 << endl;
		//		cout << "hitsFrom2in2 : " << hitsFrom2in2 << endl;
		//		cout << "hitsFrom1in2 : " << hitsFrom1in2 << endl;
		//		cout << "hitsFrom2in1 : " << hitsFrom2in1 << endl;
		//		cout << "contamination1 : " << contamination1 << endl;
		//		cout << "contamination2 : " << contamination2 << endl;
		//		cout << "purity1 : " << purity1 << endl;
		//		cout << "purity2 : " << purity2 << endl;
		//		cout << "compensation1 : " << compensation1 << endl;
		//		cout << "compensation2 : " << compensation2 << endl;


		//		cc->Update();
		//		cc->WaitPrimitive();



		//
		//		if( nbOfPeaksx != 2 ) {
		//			cout << "couldn't find two peaks !" << endl;
		//		}


				delete xNew1;
				delete xNew2;
				delete yNew1;
				delete yNew2;
				delete x;
				delete y;
				delete xNew;
				delete yNew;
				delete newXHist;
				delete s1;


			}
		}

		else if( estimatorModeArg.getValue() == "cylinder" ) {

			double hitsFrom1Side1OutsideCylinder = 0;
			double hitsFrom2Side1OutsideCylinder = 0;
			double hitsFrom2Side2OutsideCylinder = 0;
			double hitsFrom1Side2OutsideCylinder = 0;

			ThreeVector cog1 = GetCenterOfGravity(hitColType1);
			double rms1 = GetRMSDispersion(hitColType1);
			ThreeVector center1(cog1.x(),cog1.y(),nbOfPadsXYZ.at(2)/2);
			ThreeVector direc1(0,0,1);
			Cylinder *cylinder1 = new Cylinder(center1,direc1,nbOfPadsXYZ.at(2),rms1/2.0);


			ThreeVector cog2 = GetCenterOfGravity(hitColType2);
			double rms2 = GetRMSDispersion(hitColType2);
			ThreeVector center2(cog2.x(),cog2.y(),nbOfPadsXYZ.at(2)/2);
			ThreeVector direc2(0,0,1);
			Cylinder *cylinder2 = new Cylinder(center2,direc2,nbOfPadsXYZ.at(2),rms2/2.0);

			double planeXPosition = ( cog1.x()/rms1 + cog2.x()/rms2 ) / ( 1/rms1 + 1/rms2 );

			for( unsigned int i=0 ; i<hitColType1->size() ; i++ ) {

				IntVec ijk = hitColType1->at(i)->GetIJK();
				ThreeVector ijkPos( ijk.at(0) , ijk.at(1) , ijk.at(2) );


				if( cylinder1->Contains(ijkPos) ) {

//					if( cylinder2->Contains(ijkPos) ) {
//						commonHits++;
//					}
					hitsFrom1in1++;
				}
				if( cylinder2->Contains(ijkPos) ) {
					hitsFrom1in2++;
				}
				if( !cylinder1->Contains(ijkPos) && !cylinder2->Contains(ijkPos) ) {
	//				misMatch1++;
					if( ijk.at(0) > planeXPosition ) hitsFrom1Side1OutsideCylinder ++;
					else hitsFrom2Side1OutsideCylinder ++;
				//if( planeXPosition - cog1.x() < 0 ) cout << "hits 1 on the right" << endl;
				//else  cout << "hits 1 on the left" << endl;


				}

			}
	
			for( unsigned int i=0 ; i<hitColType2->size() ; i++ ) {

				IntVec ijk = hitColType2->at(i)->GetIJK();
				ThreeVector ijkPos( ijk.at(0) , ijk.at(1) , ijk.at(2) );

				if( cylinder2->Contains(ijkPos) ) {

//					if( cylinder1->Contains(ijkPos) ) {
//						commonHits++;
//					}
					hitsFrom2in2++;
				}
				if( cylinder1->Contains(ijkPos) ) {
					hitsFrom2in1++;
				}
				if( !cylinder2->Contains(ijkPos) && !cylinder1->Contains(ijkPos) ) {
	//				misMatch2++;
					if( ijk.at(0) < planeXPosition ) hitsFrom2Side2OutsideCylinder ++;
					else hitsFrom1Side2OutsideCylinder ++;
				}

			}

			contamination1 = ( hitsFrom2in1 + hitsFrom2Side1OutsideCylinder ) / ( hitsFrom1in1 + hitsFrom2in1 + hitsFrom1Side1OutsideCylinder + hitsFrom2Side1OutsideCylinder);
			contamination2 = ( hitsFrom1in2 + hitsFrom1Side2OutsideCylinder ) / ( hitsFrom1in2 + hitsFrom2in2 + hitsFrom2Side2OutsideCylinder + hitsFrom1Side2OutsideCylinder);

			purity1 = ( hitsFrom1in1 + hitsFrom1Side1OutsideCylinder ) / ( hitColType1->size() );
			purity2 = ( hitsFrom2in2 + hitsFrom2Side2OutsideCylinder ) / ( hitColType2->size() );

			compensation1 = ( hitsFrom2in1 + hitsFrom2Side1OutsideCylinder ) / hitColType1->size();
			compensation2 = ( hitsFrom1in2 + hitsFrom1Side2OutsideCylinder ) / hitColType2->size();

			delete cylinder1;
			delete cylinder2;

		}

		analysisManager->Set("EstimatorVariables","purity1",purity1);
		analysisManager->Set("EstimatorVariables","purity2",purity2);
		analysisManager->Set("EstimatorVariables","contamination1",contamination1);
		analysisManager->Set("EstimatorVariables","contamination2",contamination2);
		analysisManager->Set("EstimatorVariables","compensation1",compensation1);
		analysisManager->Set("EstimatorVariables","compensation2",compensation2);
		analysisManager->Set("EstimatorVariables","hitsFrom1in1",hitsFrom1in1);
		analysisManager->Set("EstimatorVariables","hitsFrom2in2",hitsFrom2in2);
		analysisManager->Set("EstimatorVariables","hitsFrom2in1",hitsFrom2in1);
		analysisManager->Set("EstimatorVariables","hitsFrom1in2",hitsFrom1in2);
		analysisManager->Set("EstimatorVariables","showersFound",showersFound);


		analysisManager->Fill("EstimatorVariables");


		delete hitColType1;
		delete hitColType2;

		hitManager->EndOfEvent();

	}


	/*******************
	 * Finalize the task
	 *******************/

	lcReader->close();
	analysisManager->End();
	AnalysisManager::Kill();
	SdhcalConfig::Kill();
	delete lcReader;


	return 0;
}

