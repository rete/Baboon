/*
 *
 * HoughTransformAlgorithm.hh header template generated by fclass
 * Creation date : Thu Mar 14 21:55:13 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#ifndef HOUGHTRANSFORMALGORITHM_HH
#define HOUGHTRANSFORMALGORITHM_HH

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <map>
#include <vector>
#include <utility>
#include <algorithm>


// sdhcal includes
#include "Algorithm/AbstractAlgorithm.hh"
#include "Objects/Cluster.hh"
#include "Objects/Hit.hh"
#include "Managers/ClusteringManager.hh"
#include "Managers/HitManager.hh"
#include "Managers/TrackManager.hh"
#include "Reconstruction/Tag.hh"
#include "Reconstruction/Builders/TrackCollectionBuilder.hh"
#include "Utilities/ReturnValues.hh"

namespace baboon {




	/*!
	 * Class HoughTransformAlgorithm.
	 * Inherit from 'AbstractAlgorithm' base class.
	 */

	class HoughTransformAlgorithm : public AbstractAlgorithm {


		/*!
		 * Hough tag used in HoughTransformAlgorithm to tag candidate clusters
		 */

		enum HoughTag {
			fGood,
			fBad
		};


		/*!
		 *
		 * @brief Cluster structure for Hough Transform algorithm. Allows to tag the clusters while processing
		 *
		 */
		typedef struct {

			public :
				Cluster *cluster;
				std::vector<int> rhox;
				std::vector<int> rhoy;
				HoughTag tagx;
				HoughTag tagy;
				HoughTag finalTag;

		} HoughCluster;


		typedef std::vector<HoughCluster*> HoughClusterCollection;


		public :
			/*!
			 *
			 * @brief Default Constructor
			 *
			 */
			HoughTransformAlgorithm();
			/*!
			 *
			 * @brief Default Destructor
			 *
			 */
			virtual ~HoughTransformAlgorithm();


		protected :

			/*!
			 *
			 * @brief Initialize the algorithm, i.e by initializing specific variables
			 *
			 */
			virtual void Init();


			/*!
			 *
			 * @brief Execute the algorithm
			 *
			 */
			virtual void Execute();


			/*!
			 *
			 * @brief Finalize the algorithm
			 *
			 */
			virtual void End();


			/*!
			 *
			 * @brief Allow to check if everything is well set in the algorithm before starting it
			 *
			 */
			virtual Return CheckConsistency();


			/*!
			 *
			 * @brief Delete the Hough parameter space ( 2D array )
			 *
			 */
			void DeleteHoughSpace();


			/*!
			 *
			 * @brief Allocate the Hough parameter space ( 2D array )
			 *
			 */
			void AllocateHoughSpace();


			// Algorithm parameters
			int thetaMax;
			int rMax;
			int clusterSizeLimit;
			int minimumBinning;
			int deltaPosMax;
			int trackSegmentMinimumSize;
			int maximumDistanceBetweenHitsInPlane;
			int maximumDistanceBetweenHitsForLayers;

			int ** houghSpaceX;
			int ** houghSpaceY;

	};

}

#endif  // HOUGHTRANSFORMALGORITHM_HH
