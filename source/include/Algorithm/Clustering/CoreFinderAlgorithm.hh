/*
 *
 * CoreFinderAlgorithm.hh header template generated by fclass
 * Creation date : Mon Apr  8 17:46:03 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#ifndef COREFINDERALGORITHM_HH
#define COREFINDERALGORITHM_HH

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "Algorithm/AbstractAlgorithm.hh"
#include "Managers/ClusteringManager.hh"
#include "Exception.hh"
#include "Managers/HitManager.hh"
#include "Objects/Hit.hh"


namespace sdhcal {


	/*!
	 * Class CoreFinderAlgorithm.
	 * Inherit from 'AbstractAlgorithm' base class. \n
	 * Find the core hits and tag them as fCore.
	 */

	class CoreFinderAlgorithm : public AbstractAlgorithm {

		protected :

			/*! HitCollection needed for the algorithm */
			HitCollection *hitCollection;

			/*! Threshold accumulation. cf. Execute() */
			std::vector<int> threshCountVec;

			/*! Threshold to be reached to tag a hit as core hit. cf. Execute() */
			int coreCountThreshold;

			/*! Distance to cluster the core hit. cf. Execute() */
			int distance;

			double minimumThresholdConcentration;

		public :

			/*! Default Constructor */
			CoreFinderAlgorithm();

			/*! Default Destructor */
			virtual ~CoreFinderAlgorithm();

			/*! Initialize the algorithm, i.e by initializing specific variables */
			virtual void Init();

			/*! Execute the algorithm */
			virtual void Execute() ;

			/*! Finalize the algorithm*/
			virtual void End();

			/*! Allow to check if everything is well set in the algorithm before starting it */
			virtual Return CheckConsistency();

			/*! Set the hit collection needed for the algorithm */
			inline void SetHitCollection( HitCollection *hitCol )
				{ hitCollection = hitCol; }

			/*! Set the distance between two pads to clusterize core hits */
			inline void SetDistance( int d )
				{ distance = d; }


	};

}

#endif  // COREFINDERALGORITHM_HH