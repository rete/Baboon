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

// sdhcal includes
#include "Geometry/ThreeVector.hh"
#include "Objects/Hit.hh"
#include "Config/SdhcalConfig.hh"




namespace baboon {


class Hit;
typedef std::vector<Hit*> HitCollection;

	/*!
	 * Class SDHCALCluster.
	 */

	class Cluster {


	protected :

		ThreeVector position;

		std::string codingPattern;

		HitCollection *hitCollection;

		std::vector<Cluster*> *associatedClusters;

		Tag fClusterTag;

		void ComputePosition();

	public :
		/*! Default Constructor */
		Cluster();

		/*! Default Destructor */
		virtual ~Cluster();

		void SetHitCollection(HitCollection *hitCol);

		const ThreeVector& GetPosition() {
			this->ComputePosition();
			return position;
		}

		inline HitCollection *GetHitCollection() const
			{ return hitCollection; }

		bool IsIsolatedFromClusters(const std::vector<Cluster*>* clusters);

		inline int GetClusterSize()
			{ return hitCollection->size(); }

		void SetClusterTagRecursive(Tag clustTag);

		inline void SetClusterTag(Tag clustTag)
			{ fClusterTag = clustTag; }

		inline Tag GetClusterTag() const
			{ return fClusterTag; }

		void AssociateCluster( Cluster* );

		std::vector<Cluster*> *GetAssociatedClusters()
			{ return associatedClusters; }

		bool IsAssociatedToCluster( Cluster* );

		bool ContainsHit( Hit *hit );

		void AddHit( Hit *hit );

//		void AddHitRecursive( Hit *hit );

		void RemoveHit( Hit *hit );

		void MergeClusters( Cluster *cl );


	};

	typedef std::vector<Cluster*> ClusterCollection;


}  //  end namespace sdhcal


#endif  // CLUSTER_HH
