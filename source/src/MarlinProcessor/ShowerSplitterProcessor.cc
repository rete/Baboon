  /// \file ShowerSplitterProcessor.cc
/*
 *
 * ShowerSplitterProcessor.cc source template generated by fclass
 * Creation date : jeu. mai 9 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "MarlinProcessor/ShowerSplitterProcessor.hh"


ShowerSplitterProcessor aShowerSplitterProcessor;


using namespace std;
using namespace baboon;

ShowerSplitterProcessor::ShowerSplitterProcessor()
	: marlin::Processor("ShowerSplitterProcessor") {

	  _description = "ShowerSplitterProcessor for shower splitting in SDHCAL";

	  // register steering parameters: name, description, class-variable, default value

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
				     string("ShowerSplitter.root"));

	  registerProcessorParameter("decoderString" ,
				     "decoder string for cell ID decoder" ,
				     decoderString,
				     string("M:3,S-1:3,I:9,J:9,K-1:6"));

	  registerProcessorParameter("collectionName" ,
				     "collection name for SDHCAL hits" ,
				     collectionName,
				     string("HCALBarrel"));

}

ShowerSplitterProcessor::~ShowerSplitterProcessor() {}


void ShowerSplitterProcessor::init() {

	/************************
	 * SDHCAL configurations
	 ************************/

	SdhcalConfig *config = SdhcalConfig::GetInstance();
	cout << "Loading configuration file : " << configFileName << endl;
	config->LoadFile( configFileName );

	/***********************************************************
	 * Define all algorithms, add them in the algorithm manager
	 ***********************************************************/

	AlgorithmManager *algorithmManager = AlgorithmManager::GetInstance();

	algorithmManager->SetConfigFileName( algorithmFileName );

	// Add the Hough Transform Algorithm for track reconstruction within the sdhcal
	algorithmManager->RegisterAlgorithm( new HoughTransformAlgorithm() );

	// Add isolation tagging algorithm
	algorithmManager->RegisterAlgorithm( new IsolationTaggingAlgorithm() );

	// Add clustering (2D) algorithm
	algorithmManager->RegisterAlgorithm( new ClusteringAlgorithm() );

	// Add core finder algorithm
	algorithmManager->RegisterAlgorithm( new CoreFinderAlgorithm() );

	// Add cone beginning algorithm
	algorithmManager->RegisterAlgorithm( new ConeBeginningAlgorithm() );

	// Add pca algorithm
//	algorithmManager->RegisterAlgorithm( new PrincipalComponentAnalysis() );

	cout << "Loading configuration file for algorithm : " << algorithmFileName << endl;
	algorithmManager->Initialize();

	/****************************************************************
	 * Define the analysis manager and give it the root output file
	 ****************************************************************/

	AnalysisManager *analysisManager = AnalysisManager::GetInstance();
	cout << "ROOT output file : " << rootOutputFile << endl;
	analysisManager->SetRootFileName( rootOutputFile );
	analysisManager->Init();

}

void ShowerSplitterProcessor::processRunHeader( LCRunHeader* run ) {

}

void ShowerSplitterProcessor::processEvent( LCEvent * evt ) {

	baboon::CaloHitCreator *caloHitCreator = new baboon::CaloHitCreator();
	caloHitCreator->SetDecoderString( decoderString );
	caloHitCreator->SetCollectionName( collectionName );
	caloHitCreator->CreateCaloHits( evt );

	HitManager *hitManager = HitManager::GetInstance();
	hitManager->BuildVolumeMap();

	AnalysisManager *analysisManager = AnalysisManager::GetInstance();

	// just to keep tree compatibility with the old version. Will be removed in the future
	int evtNb = evt->getEventNumber();
	cout << "event " << evtNb <<  endl;
	analysisManager->Set("SplitterVariables","event",evtNb);
	analysisManager->Set("SplitterVariables","nbOfClusters",0);
	analysisManager->Set("SplitterVariables", "nbOfTracks" , 0 );

	vector<int> I;
	vector<int> J;
	vector<int> K;

	vector<int> ITrack;
	vector<int> JTrack;
	vector<int> KTrack;

	vector<int> IThrust;
	vector<int> JThrust;
	vector<int> KThrust;

	vector<int> ICore;
	vector<int> JCore;
	vector<int> KCore;

	vector<int> ITrackExtr;
	vector<int> JTrackExtr;
	vector<int> KTrackExtr;

	vector<double> ITrackFitFirst;
	vector<double> JTrackFitFirst;
	vector<double> KTrackFitFirst;

	vector<double> ITrackFitSecond;
	vector<double> JTrackFitSecond;
	vector<double> KTrackFitSecond;

	vector<int> IIsolated;
	vector<int> JIsolated;
	vector<int> KIsolated;

	vector<double> Chi2Vec;
	vector<double> isolationWeights;

	AlgorithmManager *algorithmManager = AlgorithmManager::GetInstance();
	ClusteringManager *clusteringManager = ClusteringManager::GetInstance();
	TrackManager *trackManager = TrackManager::GetInstance();
	CoreManager *coreManager = CoreManager::GetInstance();

	HitCollection * hitCollection = hitManager->GetHitCollection();
	ClusterCollection *clusterCollection = clusteringManager->GetCluster2D();

//	cout << "nb of clusters : " << clusterCollection->size() << endl;
//	cout << "nb of hits : " << hitCollection->size() << endl;



	if( algorithmManager->AlgorithmIsRegistered("IsolationTaggingAlgorithm") ) {

		cout << "IsolationTaggingAlgorithm found" << endl;
		IsolationTaggingAlgorithm* isolationAlgo = (IsolationTaggingAlgorithm *) algorithmManager->GetAlgorithm("IsolationTaggingAlgorithm");
		isolationAlgo->Process();
		isolationWeights = isolationAlgo->GetIsolationWeights();
	}

	if( algorithmManager->AlgorithmIsRegistered("CoreFinderAlgorithm") ) {

		cout << "CoreFinderAlgorithm found" << endl;
		CoreFinderAlgorithm *coreFinder = (CoreFinderAlgorithm *) algorithmManager->GetAlgorithm("CoreFinderAlgorithm");
		coreFinder->Process();
	}

	if( algorithmManager->AlgorithmIsRegistered("ClusteringAlgorithm") ) {

		cout << "ClusteringAlgorithm found" << endl;
		ClusterCollection *clustCol = new ClusterCollection();
		ClusteringAlgorithm* clusteringAlgo = (ClusteringAlgorithm *) algorithmManager->GetAlgorithm("ClusteringAlgorithm");
		clusteringAlgo->SetClusteringMode( fClustering3D );
		clusteringAlgo->SetTaggingMode( fClusterTagMode );
		clusteringAlgo->AddHitTagToCluster( fCore );
		clusteringAlgo->SetClusterCollection( clustCol );
		clusteringAlgo->Process();

		for( unsigned int i=0 ; i<clustCol->size() ; i++ ) {

			HitCollection *hitCol = clustCol->at(i)->GetHitCollection();
			Core *core = new Core();
			for( unsigned int j=0 ; j<hitCol->size() ; j++ ) {
				core->AddHit( hitCol->at(j) );
			}
			Return ret = CoreManager::GetInstance()->AddCore( core );
			if( !ret.OK )
				cout << ret.message << endl;
		}
		clustCol->clear();
		delete clustCol;
	}

	if( algorithmManager->AlgorithmIsRegistered("ConeBeginningAlgorithm") ) {

		cout << "ConeBeginningAlgorithm found" << endl;
		ConeBeginningAlgorithm *coneBeginning = (ConeBeginningAlgorithm *) algorithmManager->GetAlgorithm("ConeBeginningAlgorithm");
		coneBeginning->Process();
	}

	if( algorithmManager->AlgorithmIsRegistered("ClusteringAlgorithm") ) {

		cout << "ClusteringAlgorithm found" << endl;
		ClusterCollection *clustCol = new ClusterCollection();
		ClusteringAlgorithm* clusteringAlgo = (ClusteringAlgorithm *) algorithmManager->GetAlgorithm("ClusteringAlgorithm");
		clusteringAlgo->SetClusteringMode( fClustering2D );
		clusteringAlgo->SetTaggingMode( fClusterTagMode );
		clusteringAlgo->SetClusterCollection( clustCol );
		clusteringAlgo->Process();
		for( unsigned int i=0 ; i<clustCol->size() ; i++ ) {
			Return ret = ClusteringManager::GetInstance()->AddCluster( clustCol->at(i) );
			if( !ret.OK )
				cout << ret.message << endl;
		}
		clustCol->clear();
		delete clustCol;
	}

	if( algorithmManager->AlgorithmIsRegistered("HoughTransformAlgorithm") ) {

		cout << "HoughTransformAlgorithm found" << endl;
		HoughTransformAlgorithm *houghTransform = (HoughTransformAlgorithm*) algorithmManager->GetAlgorithm("HoughTransformAlgorithm");
		houghTransform->Process();
	}




	for(unsigned int j=0 ; j<hitCollection->size() ; j++) {

		IntVec ijk = hitCollection->at(j)->GetIJK();

		if( hitCollection->at(j)->GetHitTag() == fTrack ) {
			ITrack.push_back( ijk.at(0) );
			JTrack.push_back( ijk.at(1) );
			KTrack.push_back( ijk.at(2) );
		}
		else if (hitCollection->at(j)->GetHitTag() == fTrackExtremity) {
			ITrackExtr.push_back( ijk.at(0) );
			JTrackExtr.push_back( ijk.at(1) );
			KTrackExtr.push_back( ijk.at(2) );
		}
		else if( hitCollection->at(j)->GetHitTag() == fCore ) {
			ICore.push_back( ijk.at(0) );
			JCore.push_back( ijk.at(1) );
			KCore.push_back( ijk.at(2) );
		}
		else if( hitCollection->at(j)->GetHitTag() == fIsolated ) {
			IIsolated.push_back( ijk.at(0) );
			JIsolated.push_back( ijk.at(1) );
			KIsolated.push_back( ijk.at(2) );
		}
		else {
			I.push_back( ijk.at(0) );
			J.push_back( ijk.at(1) );
			K.push_back( ijk.at(2) );
		}
	}

	analysisManager->Set("SplitterVariables","Hitx",&I);
	analysisManager->Set("SplitterVariables","Hity",&J);
	analysisManager->Set("SplitterVariables","Hitz",&K);

	analysisManager->Set("SplitterVariables","HitxTrack",&ITrack);
	analysisManager->Set("SplitterVariables","HityTrack",&JTrack);
	analysisManager->Set("SplitterVariables","HitzTrack",&KTrack);

	analysisManager->Set("SplitterVariables","HitxThrust",&IThrust);
	analysisManager->Set("SplitterVariables","HityThrust",&JThrust);
	analysisManager->Set("SplitterVariables","HitzThrust",&KThrust);

	analysisManager->Set("SplitterVariables","HitxCore",&ICore);
	analysisManager->Set("SplitterVariables","HityCore",&JCore);
	analysisManager->Set("SplitterVariables","HitzCore",&KCore);

	analysisManager->Set("SplitterVariables","HitxTrackExtr",&ITrackExtr);
	analysisManager->Set("SplitterVariables","HityTrackExtr",&JTrackExtr);
	analysisManager->Set("SplitterVariables","HitzTrackExtr",&KTrackExtr);

	analysisManager->Set("SplitterVariables","ITrackFitFirst",&ITrackFitFirst);
	analysisManager->Set("SplitterVariables","JTrackFitFirst",&JTrackFitFirst);
	analysisManager->Set("SplitterVariables","KTrackFitFirst",&KTrackFitFirst);

	analysisManager->Set("SplitterVariables","ITrackFitSecond",&ITrackFitSecond);
	analysisManager->Set("SplitterVariables","JTrackFitSecond",&JTrackFitSecond);
	analysisManager->Set("SplitterVariables","KTrackFitSecond",&KTrackFitSecond);

	analysisManager->Set("SplitterVariables","IIsolated",&IIsolated);
	analysisManager->Set("SplitterVariables","JIsolated",&JIsolated);
	analysisManager->Set("SplitterVariables","KIsolated",&KIsolated);

	analysisManager->Set("SplitterVariables","Chi2Vec",&Chi2Vec);
	analysisManager->Set("SplitterVariables","isolationWeights",&isolationWeights);

	analysisManager->Fill("SplitterVariables");

	clusteringManager->ClearAllContent();
	hitManager->ClearAllContent();
	coreManager->ClearAllContent();
	trackManager->ClearAllContent();

	delete caloHitCreator;
}


void ShowerSplitterProcessor::check( LCEvent *evt ) {

}

void ShowerSplitterProcessor::end() {

	AnalysisManager::GetInstance()->End();

	AnalysisManager::Kill();
	AlgorithmManager::Kill();
	HitManager::Kill();
	ClusteringManager::Kill();
	CoreManager::Kill();
	SdhcalConfig::Kill();

}


