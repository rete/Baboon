  /// \file IsolatedHitMergingAlgorithm.hh
/*
 *
 * IsolatedHitMergingAlgorithm.hh header template generated by fclass
 * Creation date : mer. d�c. 11 2013
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
 * @author : Remi Ete
 * @version
 * @copyright
 *
 *
 */


#ifndef ISOLATEDHITMERGINGALGORITHM_HH
#define ISOLATEDHITMERGINGALGORITHM_HH

#include "Algorithm/AbstractAlgorithm.hh"

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "Detector/Calorimeter.hh"
#include "Managers/ClusteringManager.hh"
#include "Objects/Cluster.hh"

namespace baboon {

	/**
	 * @brief  IsolatedHitMergingAlgorithm class
	 * Inherits from base class AbstractAlgorithm
	 */
	class IsolatedHitMergingAlgorithm : public AbstractAlgorithm {

	public:

		/**
		 * @brief  Default constructor
		 */
		IsolatedHitMergingAlgorithm();

		/**
		 * @brief  Default destructor
		 */
		virtual ~IsolatedHitMergingAlgorithm();

		/**
		 *
		 */
		void SetCalorimeter( Calorimeter *calo )
			{ calorimeter = calo; }

	protected:

		/**
		 * @brief Initialize the algorithm, i.e by initializing specific variables
		 */
		Return Init();

		/**
		 * @brief Allow to check if everything is well set in the algorithm before starting it
		 */
		Return CheckConsistency();

		/**
		 * @brief Execute the algorithm
		 */
		Return Execute();

		/**
		 * @brief Finalize the algorithm
		 */
		Return End();

		/**
		 *
		 */
		double DistanceToCluster( const ThreeVector &pos , Cluster *cluster );



		Calorimeter *calorimeter;

	};  // class

}  // namespace 

#endif  //  ISOLATEDHITMERGINGALGORITHM_HH
