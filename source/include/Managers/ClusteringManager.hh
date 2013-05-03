/*
 *
 * ClusteringManager.hh header template generated by fclass
 * Creation date : Fri Mar 15 18:06:44 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#ifndef CLUSTERINGMANAGER_HH
#define CLUSTERINGMANAGER_HH

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>


#include "Objects/Cluster.hh"
#include "Objects/Hit.hh"
#include "Managers/HitManager.hh"

namespace sdhcal {

	/*!
	 * Class ClusteringManager.
	 */

	class ClusteringManager {

	protected :

		/*! Default Constructor */
		ClusteringManager();

		/*! Default Destructor */
		virtual ~ClusteringManager();

		static ClusteringManager *instance;

		std::string codingPattern;

		ClusterCollection *clusters3D;

		ClusterCollection *clusters2D;


	public :

		static ClusteringManager *GetInstance();

		static void Kill();

		ClusterCollection *GetCluster3D();

		ClusterCollection *GetCluster2D();


	};

}

#endif  // CLUSTERINGMANAGER_HH
