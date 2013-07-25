  /// \file BaboonMonitoring.hh
/*
 *
 * BaboonMonitoring.hh header template generated by fclass
 * Creation date : jeu. juil. 11 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef BABOONMONITORING_HH
#define BABOONMONITORING_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector>
#include <map>

#include "TApplication.h"
#include <TGeometry.h>
#include <TGeoMaterial.h>
#include <TGeoManager.h>
#include <TEveGeoNode.h>

#include "Xml/tinyxml.h"
#include "Xml/XmlHelper.h"

#include "Reconstruction/Tag.hh"
#include "Objects/HitCollection.hh"
#include "Managers/ShowerManager.hh"
#include "Managers/DetectorManager.hh"
#include "Detector/Calorimeter.hh"

namespace baboon {

//	/*!
//	 * @brief Display mode enum
//	 */
//	enum HitDisplayMode {
//
//		kDisplayThresholds,
//		kDisplayShowers,
//		kDisplayTags,
//		kDisplayUniform
//	};


	/*
	* Class BaboonMonitoring
	*/

	class BaboonMonitoring {

		public:

			/*!
			 *
			 *
			 *
			 */
			static BaboonMonitoring *GetInstance();

			/*!
			 *
			 *
			 *
			 */
			static void Kill();

			/*!
			 *
			 *
			 *
			 */
			void Pause() const;

			/*!
			 *
			 *
			 *
			 */
			void ViewEvent();


		protected:

			/*!
			 *
			 * Default Constructor
			 *
			 * */
			BaboonMonitoring();

			/*!
			 *
			 *  Default Destructor
			 *
			 */
			~BaboonMonitoring();

			/*!
			 *
			 *
			 *
			 */
			void InitializeEve(Char_t transparency = 70);

			/*!
			 *
			 *
			 *
			 */
			void InitializeDetectors(TGeoVolume *pMainDetectorVolume, TGeoMedium *pSubDetectorMedium, Char_t transparency);

			/*!
			 *
			 *
			 *
			 */
//			TEveElement *VisualizeHitCollection(const HitCollection *hitCollection, TEveElement *parent, const DisplayMode displayMode );

			/*!
			 *
			 *
			 *
			 */
			int GetHitColor( Hit *hit , const HitDisplayMode displayMode );

			static BaboonMonitoring *instance;

			TApplication *rootApplication;
			static std::map<Tag,int> tagToColorMap;

			static bool openEveEvent;
			static int eventDisplayCounter;
			static bool isTEveInitialized;

	};  // class

}  // namespace 

#endif  //  BABOONMONITORING_HH
