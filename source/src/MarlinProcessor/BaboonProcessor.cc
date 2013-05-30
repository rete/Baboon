  /// \file BaboonProcessor.cc
/*
 *
 * BaboonProcessor.cc source template generated by fclass
 * Creation date : lun. mai 27 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "MarlinProcessor/BaboonProcessor.hh"

using namespace std;
using namespace EVENT;

namespace baboon {

	BaboonProcessor::BaboonProcessor( const std::string &processorName )
		: marlin::Processor( processorName ) {



		  registerProcessorParameter( "BABOON_HOME" ,
					      "Path to Baboon directory install",
					      baboonHome,
					      string("/home/remi/ilcsoft/SDHCAL/Baboon") );

		  string cfgFileName = string("/home/remi/ilcsoft/SDHCAL/Baboon/config/SDHCAL.cfg");
		  registerProcessorParameter("SDHCAL_cfg" ,
					     "SDHCAL configurations" ,
					     configFileName,
					     cfgFileName);

		  string algoFileName = string("/home/remi/ilcsoft/SDHCAL/Baboon/config/Algorithm.cfg");
		  registerProcessorParameter("Algorithm_cfg" ,
					     "Algorithm configurations" ,
					     algorithmFileName,
					     algoFileName);

		  registerProcessorParameter("rootOutputFile" ,
					     "root outputfile" ,
					     rootOutputFile,
					     string(""));

		  registerProcessorParameter("decoderString" ,
					     "decoder string for cell ID decoder" ,
					     decoderString,
					     string("M:3,S-1:3,I:9,J:9,K-1:6"));

		  registerProcessorParameter("collectionName" ,
					     "collection name for SDHCAL hits" ,
					     collectionName,
					     string("HCALBarrel"));

		  registerProcessorParameter("graphicalEnvironment" ,
					     "Set the GUI ON or OFF while processing" ,
					     graphicalEnvironment,
					     false);

	}

	BaboonProcessor::~BaboonProcessor() {

		if( graphicalEnvironment ) {
			gApplication->Terminate();
			delete gApplication;
		}
	}





	void BaboonProcessor::init() {

		algorithmManager = AlgorithmManager::GetInstance();

		if( graphicalEnvironment ) {

			string appName("BaboonProcessorGUI");
			if( gApplication == 0 )
				gApplication = new TRint(appName.c_str(),0,0);
		}


		/*
		 *
		 * Function implemented by the user. Called after
		 * loading the graphical environment if there's one.
		 * Called before loading the managers and load the graphics.
		 *
		 */
		BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , this->Init() );


		this->LoadManagers();

		if( graphicalEnvironment ) {

			calorimeter = SDHCALPrototype::GetInstance();
			calorimeter->BuildGeometry();
			TGeoVolume *caloVolume = SDHCALPrototype::GetGeoManager()->GetTopVolume();
			caloVolume->Draw("ogl");
			TGLViewer *viewer = (TGLViewer*) gPad->GetViewer3D("ogl");
			viewer->SetClearColor( TColor::GetColor(133,194,163) );
			gPad->GetCanvas()->SetWindowSize(100,100);
		}

	}



	void BaboonProcessor::processRunHeader( LCRunHeader* run ) {

		BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , this->ProcessRunHeader( run ) );
	}


	void BaboonProcessor::processEvent( EVENT::LCEvent * evt ) {

		// Load the event in the baboon framework and initialize the objects.
		this->LoadEvent( evt );

		// Load the hits on the current scene id the graphical environment is set.
		if( graphicalEnvironment )
			calorimeter->LoadHitCollection( hitManager->GetHitCollection() );


		/*
		 *
		 * Process an event in the BABOON framework.
		 * Must be implemented by the user.
		 *
		 */
		BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , this->ProcessEvent( evt->getEventNumber() ) );


		if( graphicalEnvironment ) {

			gPad->Update();
			gPad->WaitPrimitive();
		}


		this->ClearAllContent();
	}



	void BaboonProcessor::check( LCEvent * evt ) {

		BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , this->Check( evt ) );
	}



	void BaboonProcessor::end() {

		BABOON_THROW_RESULT_IF( BABOON_SUCCESS() , != , this->End() );

		AnalysisManager::GetInstance()->End();
		AnalysisManager::Kill();
		AlgorithmManager::Kill();
		HitManager::Kill();
		ClusteringManager::Kill();
		CoreManager::Kill();
		ShowerManager::Kill();
		SdhcalConfig::Kill();
		SDHCALPrototype::Kill();
	}



	void BaboonProcessor::LoadManagers() {

		/************************
		 * SDHCAL configurations
		 ************************/

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

		if( !rootOutputFile.empty() )
			cout << "ROOT output file : " << rootOutputFile << endl;
		analysisManager->SetRootFileName( rootOutputFile );
		analysisManager->Init();


		/**************************************************************
		 * Load the algorithm manager with the algorithm configuration
		 **************************************************************/

		algorithmManager->SetConfigFileName( algorithmFileName );
		algorithmManager->Initialize();
		cout << "end of load managers!" << endl;

	}




	void BaboonProcessor::LoadEvent( LCEvent *evt ) {


		caloHitCreator = new CaloHitCreator();
		caloHitCreator->SetDecoderString( decoderString );
		caloHitCreator->SetCollectionName( collectionName );
		caloHitCreator->CreateCaloHits( evt );

		HitManager *hitManager = HitManager::GetInstance();
		hitManager->BuildVolumeMap();

	}




	void BaboonProcessor::ClearAllContent() {

		clusteringManager->ClearAllContent();
		hitManager->ClearAllContent();
		coreManager->ClearAllContent();
		trackManager->ClearAllContent();
		showerManager->ClearAllContent();

		if( graphicalEnvironment )
			calorimeter->ClearCalorimeter();

		delete caloHitCreator;
	}


}  // namespace 

