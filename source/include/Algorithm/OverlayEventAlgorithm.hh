  /// \file OverlayEventAlgorithm.hh
/*
 *
 * OverlayEventAlgorithm.hh header template generated by fclass
 * Creation date : jeu. nov. 7 2013
 *
 * This file is part of XXX libraries.
 * 
 * XXX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * XXX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author : R�mi Et�
 * @version
 * @copyright
 *
 *
 */


#ifndef OVERLAYEVENTALGORITHM_HH
#define OVERLAYEVENTALGORITHM_HH

#include "Algorithm/AbstractAlgorithm.hh"

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "Objects/CaloHit.hh"
#include "Detector/Calorimeter.hh"

// lcio includes
#include "IMPL/ReconstructedParticleImpl.h"

namespace baboon {

	/**
	 * @brief  OverlayEventAlgorithm class
	 * Inherits from base class AbstractAlgorithm
	 */

	class OverlayEventAlgorithm : public AbstractAlgorithm {

	public:

		struct TrackInfo {

			Track *track;
			ThreeVector endPosition;
			ThreeVector beginPosition;
			ThreeVector backwardThrust;
			bool isPrimaryTrack;

		};

		/**
		 * @brief  Default constructor
		 */
		OverlayEventAlgorithm();

		/*!
		 * @brief  Default destructor
		 */
		virtual ~OverlayEventAlgorithm();

		/**
		 * @brief Initialize the algorithm, i.e by initializing specific variables
		 */
		virtual Return Init();


		/**
		 * @brief Execute the algorithm
		 */
		virtual Return Execute();


		/**
		 * @brief Finalize the algorithm
		 */
		virtual Return End();


		/**
		 * @brief Allow to check if everything is well set in the algorithm before starting it
		 */
		virtual Return CheckConsistency();

		/**
		 * @brief Returns true if the collections have been overlaid. To be checked after Process() method call.
		 */
		bool OverlayDone()
			{ return _overlayDone; }

		/**
		 *
		 */
		void SetCollection1( CaloHitCollection *caloHitCollection )
			{ _collection1 = caloHitCollection; }

		/**
		 *
		 */
		void SetCollection2( CaloHitCollection *caloHitCollection )
			{ _collection2 = caloHitCollection; }

		/**
		 *
		 */
		void SetCalorimeter1( Calorimeter *calo )
			{ _calorimeter1 = calo; }

		/**
		 *
		 */
		void SetCalorimeter2( Calorimeter *calo )
			{ _calorimeter2 = calo; }

		/**
		 *
		 */
		const CaloHitCollection *GetOverlaidCollection()
			{ return _overlaidCollection; }

		/**
		 *
		 */
		unsigned int GetNumberOfLostHits()
			{ return _lostHitCollection->size(); }

		/*
		 *
		 */
		bool GeneratesLCTracks()
			{ return _generatesLCTracks; }

		/**
		 *
		 */
		std::pair< IMPL::ReconstructedParticleImpl * , IMPL::ReconstructedParticleImpl * > &GetTrackPair()
			{ return _trackPair; }

	protected:

		/**
		 *
		 */
		void TranslateCollection( Calorimeter *calorimeter
								   , CaloHitCollection *collection
								   , CaloHitCollection *finalCollection
								   , const ThreeVector &vec );

		/**
		 *
		 */
		int ThresholdToInt( CaloHitThreshold fThr );

		/**
		 *
		 */
		void OverlayCollections( CaloHitCollection *collection1 , CaloHitCollection *collection2 );

		/**
		 *
		 */
		void FillTrackInfo( TrackInfo *trackInfo );

		/**
		 *
		 */
		void EraseTrackFromCollection( Track *track , CaloHitCollection *collection );

		void EraseTrackFromCollection( TrackInfo *trackInfo , CaloHitCollection *collection );

		ThreeVector FindShowerEnteringPoint(CaloHitCollection *collection);



		bool _overlayDone;
		bool _generatesLCTracks;
		CaloHitCollection *_collection1;
		CaloHitCollection *_collection2;
		CaloHitCollection *_collectionToOverlay1;
		CaloHitCollection *_collectionToOverlay2;
		CaloHitCollection *_lostHitCollection;
		CaloHitCollection *_overlaidCollection;
		Calorimeter *_calorimeter1;
		Calorimeter *_calorimeter2;
		std::pair< IMPL::ReconstructedParticleImpl * , IMPL::ReconstructedParticleImpl * > _trackPair;

		// Algorithm parameters
		bool _useTrackInfo;             // use the primary track info or not
		int _separationDistance;        // in pad distance
		std::string _particleType1;     // can be "charged" or "neutral"
		std::string _particleType2;
		double _inputEnergy1;
		double _inputEnergy2;

	};  // class

}  // namespace 

#endif  //  OVERLAYEVENTALGORITHM_HH
