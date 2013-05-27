  /// \file myProc.hh
/*
 *
 * myProc.hh header template generated by fclass
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


#ifndef myProc_HH
#define myProc_HH

#include "marlin/Processor.h"
#include "EVENT/CalorimeterHit.h"
#include "lcio.h"

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 


#include "Config/SdhcalConfig.hh"
#include "Algorithm/AlgorithmHeaders.hh"
#include "Managers/AlgorithmManager.hh"
#include "Managers/AnalysisManager.hh"
#include "Managers/CoreManager.hh"
#include "Managers/HitManager.hh"
#include "Managers/CoreManager.hh"
#include "Managers/ShowerManager.hh"
#include "Managers/ClusteringManager.hh"
#include "Managers/TrackManager.hh"

#include "MarlinProcessor/CaloHitCreator.hh"

/* 
 * Class myProc
 * Inherits from base class marlin::Processor
 */ 

class myProc : public marlin::Processor {

	public:

		/*!
		 *
		 * @brief Default Constructor
		 *
		 */
		myProc();

		/*!
		 *
		 * @brief Default Destructor
		 *
		 */
		virtual ~myProc();

		/*!
		 *
		 * @brief Return a new processor. Called by Marlin to execute the processor
		 *
		 */
		virtual Processor *newProcessor()
			{ return new myProc(); }

		/*!
		 *
		 * @brief Called by Marlin. Called at the begin of the job before anything is read.
		 * Use to initialize the processor, e.g. book histograms.
		 *
		 */
		virtual void init();

		/*!
		 *
		 * @brief Called by Marlin. Called for every run.
		 *
		 */
		virtual void processRunHeader( LCRunHeader* run );

		/*!
		 *
		 * @brief Called by Marlin. Called for every event - the working horse.
		 *
		 */
		virtual void processEvent( LCEvent * evt );

		/*!
		 *
		 * @brief Called by Marlin. Check the event
		 *
		 */
		virtual void check( LCEvent * evt );

		/*!
		 *
		 * @brief Called by Marlin. Called after data processing for clean up.
		 *
		 */
		virtual void end();

		/*!
		 *
		 * @brief Laod the baboon objects. Called at beginning of event in processEvent()
		 *
		 */
		void LoadEvent( LCEvent *evt );

		/*!
		 *
		 * @brief Clear all the baboon content. Called at the end of event in processEvent()
		 *
		 */
		void ClearAllContent();

		/*!
		 *
		 * @brief Load the managers. Called in init() function.
		 *
		 */
		void LoadManagers();

	protected:

		std::string baboonHome;
		std::string configFileName;
		std::string algorithmFileName;
		std::string rootOutputFile;
		std::string decoderString;
		std::string collectionName;

		baboon::CaloHitCreator *caloHitCreator;

		baboon::HitManager        *hitManager;
		baboon::ClusteringManager *clusteringManager;
		baboon::CoreManager       *coreManager;
		baboon::ShowerManager     *showerManager;
		baboon::TrackManager      *trackManager;
		baboon::AnalysisManager   *analysisManager;
		baboon::AlgorithmManager  *algorithmManager;


};  // class 


#endif  //  myProc_HH