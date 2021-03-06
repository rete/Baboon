  /// \file DetectorManager.hh
/*
 *
 * DetectorManager.hh header template generated by fclass
 * Creation date : ven. juin 28 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef DETECTORMANAGER_HH
#define DETECTORMANAGER_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 

#include "Detector/Detector.hh"
#include "Detector/Calorimeter.hh"
#include "Utilities/ReturnValues.hh"
#include "Utilities/Internal.hh"

// root includes
#include "TGeoManager.h"

// gear includes
#include "gear/GEAR.h"
#include "gear/GearMgr.h"
#include "gearxml/GearXML.h"

namespace baboon {

	/*!
	 *
	 * Class DetectorManager
	 *
	 */

	class DetectorManager {

	  public:

		/*!
		 *
		 * Return a single instance of this class. Singleton pattern
		 *
		 */
		static DetectorManager *GetInstance();

		/*!
		 *
		 * Kill the instance carefully. Singleton pattern
		 *
		 */
		static void Kill();

		/*!
		 *
		 *
		 *
		 */
		Return RegisterDetector( Detector *detector );

		/*!
		 *
		 *
		 *
		 */
		bool DetectorIsRegistered( const std::string &detectorName );

		/*!
		 *
		 *
		 *
		 */
		Return BuildGeometry( const gear::GearMgr *gearMgr );

		/*!
		 *
		 *
		 *
		 */
		Detector *GetDetector( const std::string detectorName );

		/*!
		 *
		 *
		 *
		 */
		Return Init( const gear::GearMgr *gearMgr );

		/*!
		 *
		 *
		 *
		 */
		StringVector GetDetectorList();

	  protected:

		/*!
		 *
		 * Unique instance of this class. Singleton pattern
		 *
		 */
		static DetectorManager* instance;

		/*!
		 *
		 * Protected constructor. Singleton pattern
		 *
		 */
		DetectorManager();

		/*!
		 *
		 * Protected constructor. Singleton pattern
		 *
		 */
		virtual ~DetectorManager();




	  protected:

		TGeoManager *geoManager;
		DetectorCollection *detectors;
		TGeoVolume *topVolume;
		bool geometryBuilt;


		friend class BaboonMonitoring;

	};  // class

}  // namespace 

#endif  //  DETECTORMANAGER_HH
