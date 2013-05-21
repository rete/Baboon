  /// \file GraphicalProcessor.cc
/*
 *
 * GraphicalProcessor.cc source template generated by fclass
 * Creation date : lun. mai 20 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : remi
 */


#include "GraphicalProcessor/GraphicalProcessor.hh"


using namespace std;
using namespace EVENT;

namespace baboon {


	GraphicalProcessor::GraphicalProcessor( const std::string &processorName )
		: marlin::Processor(processorName) {

	}

	GraphicalProcessor::~GraphicalProcessor() {

	}


	void GraphicalProcessor::init() {

		algorithmManager = AlgorithmManager::GetInstance();
		Return ret = this->Init();
		this->LoadManagers();
		calorimeter = SDHCALPrototype::GetInstance();
		calorimeter->BuildGeometry();
	}


	void GraphicalProcessor::processRunHeader( LCRunHeader* run ) {

		this->ProcessRunHeader( run );
	}


	void GraphicalProcessor::processEvent( LCEvent * evt ) {

		this->LoadEvent( evt );
		calorimeter->LoadHitCollection( hitManager->GetHitCollection() );
		this->ProcessEvent();
		TGeoVolume *caloVolume = calorimeter->GetCalorimeterVolume();
		caloVolume->Draw("ogl");
		TGLViewer *viewer = (TGLViewer*) gPad->GetViewer3D("ogl");
		viewer->SetClearColor(1);
		gPad->WaitPrimitive();
		this->ClearAllContent();
	}


	void GraphicalProcessor::check( LCEvent * evt ) {

		this->Check( evt );
	}

	void GraphicalProcessor::end() {

		this->End();

		AnalysisManager::GetInstance()->End();
		AnalysisManager::Kill();
		AlgorithmManager::Kill();
		HitManager::Kill();
		ClusteringManager::Kill();
		CoreManager::Kill();
		ShowerManager::Kill();
		SdhcalConfig::Kill();
	}

	void GraphicalProcessor::LoadManagers() {

		/************************
		 * SDHCAL configurations
		 ************************/
		cout << "debug" << endl;
		SdhcalConfig *config = SdhcalConfig::GetInstance();
		cout << "Loading configuration file : " << configFileName << endl;
		config->LoadFile( configFileName );


		/*****************************
		 * Get all manager instances
		 *****************************/

		clusteringManager = ClusteringManager::GetInstance();
		analysisManager = AnalysisManager::GetInstance();
		trackManager = TrackManager::GetInstance();
		coreManager = CoreManager::GetInstance();
		hitManager = HitManager::GetInstance();
		showerManager = ShowerManager::GetInstance();


		/****************************************************
		 * Load the analysis manager with the root file name
		 ****************************************************/

		cout << "ROOT output file : " << rootOutputFile << endl;
		analysisManager->SetRootFileName( rootOutputFile );
		analysisManager->Init();


		/**************************************************************
		 * Load the algorithm manager with the algorithm configuration
		 **************************************************************/

		algorithmManager->SetConfigFileName( algorithmFileName );
		algorithmManager->Initialize();
	}


	void GraphicalProcessor::LoadEvent( LCEvent *evt ) {

		caloHitCreator = new baboon::CaloHitCreator();
		caloHitCreator->SetDecoderString( decoderString );
		caloHitCreator->SetCollectionName( collectionName );
		caloHitCreator->CreateCaloHits( evt );

		HitManager *hitManager = HitManager::GetInstance();
		hitManager->BuildVolumeMap();
	}


	void GraphicalProcessor::ClearAllContent() {

		clusteringManager->ClearAllContent();
		hitManager->ClearAllContent();
		coreManager->ClearAllContent();
		trackManager->ClearAllContent();
		showerManager->ClearAllContent();
		calorimeter->ClearCalorimeter();
		delete caloHitCreator;
	}


}  // namespace 

