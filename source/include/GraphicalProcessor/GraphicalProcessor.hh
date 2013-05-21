  /// \file GraphicalProcessor.hh
/*
 *
 * GraphicalProcessor.hh header template generated by fclass
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


#ifndef GRAPHICALPROCESSOR_HH
#define GRAPHICALPROCESSOR_HH

#include "marlin/Processor.h"

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
#include "GraphicalProcessor/SDHCALPrototype.hh"


#include "TGLViewer.h"
#include "TPad.h"
//#include "TApplication.h"
//#include "TRint.h"


//extern TApplication *gApplication;


namespace baboon {

	/*!
	 * Class GraphicalProcessor
	 * Inherits from base class marlin::Processor
	 */

	class GraphicalProcessor : public marlin::Processor {

		public:
			/*! Default Constructor */
			GraphicalProcessor( const std::string & );
			/*! Default Destructor */
			virtual ~GraphicalProcessor();


			virtual Processor *newProcessor() = 0;

			/*!
			 *
			 * Called at the begin of the job before anything is read.
			 * Use to initialize the processor, e.g. book histograms.
			 *
			 */
			virtual void init();

			/*!
			 *
			 *  Called for every run.
			 *
			 */
			virtual void processRunHeader( EVENT::LCRunHeader* run );

			/*!
			 *
			 *  Called for every event - the working horse.
			 *
			 */
			virtual void processEvent( EVENT::LCEvent * evt );


			virtual void check( EVENT::LCEvent * evt );

			/*!
			 *
			 *  Called after data processing for clean up.
			 *
			 */
			virtual void end();


		protected:

			/*!
			 *
			 * @brief Laod the baboon objects. Called at beginning of event in processEvent()
			 *
			 */
			void LoadEvent( EVENT::LCEvent *evt );

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



			virtual Return Init() = 0;

			virtual Return ProcessRunHeader( EVENT::LCRunHeader* run  ) = 0;

			virtual Return ProcessEvent() = 0;

			virtual Return Check( EVENT::LCEvent * evt ) = 0;

			virtual Return End() = 0;

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

			SDHCALPrototype *calorimeter;




	};  // class

}  // namespace 

#endif  //  GRAPHICALPROCESSOR_HH
