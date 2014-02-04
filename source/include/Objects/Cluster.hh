/*
 *
 * Cluster.hh header template generated by fclass
 * Creation date : Sat Mar  2 02:06:41 2013
 * Copyright (c) 2012 CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#ifndef CLUSTER_HH
#define CLUSTER_HH

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <stdexcept>
#include <algorithm>
#include <map>

// sdhcal includes
#include "Geometry/ThreeVector.hh"
#include "Objects/HitCompositeObject.hh"
#include "Utilities/Globals.hh"
#include "Utilities/ReturnValues.hh"




namespace baboon {

	enum ClusterType {
		fCluster2D,
		fCluster3D
	};

	enum PositionComputation {
		fComputeCell,
		fComputePosition
	};

	/*!
	 * Class SDHCALCluster.
	 */

	class Cluster : public HitCompositeObject , public TypedObject {


		protected :

			ThreeVector position;
			BaseTag fClusterTag;
			ClusterType fType;

			/*!
			 *
			 *
			 *
			 */
			void ComputePosition( PositionComputation computation = fComputeCell );

		public :

			/*!
			 *
			 * Default Constructor
			 *
			 */
			Cluster();

			/*!
			 *
			 * Default Destructor
			 *
			 */
			virtual ~Cluster();

			/*!
			 *
			 *
			 *
			 */
			const ThreeVector& GetPosition( PositionComputation computation = fComputeCell ) {
				this->ComputePosition( computation );
				return position;
			}

			/*!
			 *
			 *
			 *
			 */
			void SetPosition( const ThreeVector &pos );

			/*!
			 *
			 *
			 *
			 */
			bool IsIsolatedFromClusters( const std::vector<Cluster*>* clusters );

			/*!
			 *
			 *
			 *
			 */
			void SetClusterTagRecursive( const BaseTag &clustTag );

			/*!
			 *
			 *
			 *
			 */
			inline void SetClusterTag( const BaseTag &clustTag );

			/*!
			 *
			 *
			 *
			 */
			inline BaseTag GetClusterTag() const
				{ return fClusterTag; }

			/*!
			 *
			 *
			 *
			 */
			Return SetClusterType( const ClusterType &type );

			/*!
			 *
			 *
			 *
			 */
			inline ClusterType GetClusterType() const
				{ return fType; }

			/**
			 *
			 */
			double DistanceToMe( const ThreeVector &pos );

	};

	typedef std::vector<Cluster*> ClusterCollection;
	typedef std::map<unsigned int,ClusterCollection*> OrderedClusterCollection;


}  //  end namespace sdhcal


#endif  // CLUSTER_HH
