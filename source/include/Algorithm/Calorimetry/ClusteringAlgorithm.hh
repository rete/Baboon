  /// \file ClusteringAlgorithm.hh
/*
 *
 * ClusteringAlgorithm.hh header template generated by fclass
 * Creation date : mar. mai 7 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef CLUSTERINGALGORITHM_HH
#define CLUSTERINGALGORITHM_HH

#include "Algorithm/AbstractAlgorithm.hh"

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>


#include "Objects/CaloHit.hh"
#include "Objects/Cluster.hh"
#include "Utilities/Globals.hh"

#include "Detector/Calorimeter.hh"

namespace baboon {


	enum ClusteringMode {
		fClustering2D,
		fClustering3D
	};

	enum TaggingMode {
		fAvoidTagMode,
		fClusterTagMode
	};

	/*!
	 *
	 * @brief Class ClusteringAlgorithm
	 * Inherits from base class AbstractAlgorithm
	 *
	 */

	class ClusteringAlgorithm : public AbstractAlgorithm {

		public:

			/*!
			 *
			 * @brief Default Constructor
			 *
			 */
			ClusteringAlgorithm();

			/*!
			 *
			 * @brief Default Destructor
			 *
			 */
			virtual ~ClusteringAlgorithm();



		protected:

			ClusteringMode fClusteringMode;
			TaggingMode fTaggingMode;
			TagCollection hitTagToCluster;
			TagCollection hitTagToAvoid;
			ClusterCollection *clusterCollection;
			unsigned int clusterSizeLowerLimit;
			unsigned int neighborDistance;
			CaloHitCollection treatedHits;
			Calorimeter *calorimeter;

			/*!
			 *
			 * @brief Initialize the algorithm, i.e by initializing specific variables
			 *
			 */
			virtual Return Init();


			/*!
			 *
			 * @brief Execute the algorithm
			 *
			 */
			virtual Return Execute();


			/*!
			 *
			 * @brief Finalize the algorithm
			 *
			 */
			virtual Return End();


			/*!
			 *
			 * @brief Allow to check if everything is well set in the algorithm before starting it
			 *
			 */
			virtual Return CheckConsistency();

			/*!
			 *
			 * @brief Return true if the given tag is to be avoided while clustering
			 *
			 */
			bool AvoidTag( const BaseTag &tag );

			/*!
			 *
			 * @brief Return true if the given tag has to be kept while clustering
			 *
			 */
			bool KeepTag( const BaseTag &tag );

			/*!
			 *
			 * @brief Check the tag information to continue the algorithm or not
			 *
			 */
			bool CheckTag( const BaseTag &tag );

			/*
			 *
			 * @brief Recursive method that add a hit closed to another one.
			 * This is repeated for the added hit, etc...
			 *
			 */
			Return FindCluster( CaloHit *hit , Cluster *cluster );


		public:

			/*!
			 *
			 * @brief Set the clustering mode : 2D , 3D
			 *
			 */
			Return SetClusteringMode( ClusteringMode mode );


			/*!
			 *
			 * @brief Set the tagging mode : avoid tag or cluster them
			 *
			 */
			Return SetTaggingMode( TaggingMode mode );

			/*!
			 *
			 * @brief Set the cluster collection you want to have at the end.
			 * Must be an empty one, else delete the content before processing.
			 *
			 */
			inline void SetClusterCollection( ClusterCollection *clusterCol )
				{ clusterCollection = clusterCol; }

			/*!
			 *
			 * @brief Set the hit tag to cluster hits
			 *
			 */
			Return AddHitTagToCluster( const BaseTag &tag );

			/*!
			 *
			 * @brief Set the hit tag to be avoided while clustering
			 *
			 */
			Return AddHitTagToAvoid( const BaseTag &tag );

			/*!
			 *
			 * @brief Set a lower limit for clusters.
			 *
			 */
			Return SetClusterSizeLowerLimit( unsigned int );

			/*!
			 *
			 * @brief Set the neighbor distance while clustering
			 *
			 */
			Return SetNeighborDistance( unsigned int );

			/*!
			 *
			 *
			 *
			 */
			Return SetCalorimeter( Calorimeter *calo );




	};  // class

}  // namespace 

#endif  //  CLUSTERINGALGORITHM_HH
