  /// \file EventPreparation.hh
/*
 *
 * EventPreparation.hh header template generated by fclass
 * Creation date : dim. oct. 20 2013
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


#ifndef EVENTPREPARATION_HH
#define EVENTPREPARATION_HH

#include "Algorithm/AbstractAlgorithm.hh"

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>

#include <omp.h>
#include "Managers/DetectorManager.hh"
#include "Detector/Detector.hh"
#include "Detector/Calorimeter.hh"

namespace baboon {

	/*!
	*
	* @brief  EventPreparation class
	* Inherits from base class AbstractAlgorithm
	*
	*/

	class EventPreparation : public AbstractAlgorithm {

		public:

			/*!
			*
			* @brief  Default constructor
			*
			*/
			EventPreparation();

			/*!
			*
			* @brief  Default destructor
			*
			*/
			virtual ~EventPreparation();

		protected:

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
			 *
			 *
			 */
			Return CalculateCaloHitDensities( Calorimeter *calorimeter );


			bool shouldCalculateDensities;
			int densityX;
			int densityY;
			int densityZ;


};  // class 

}  // namespace 

#endif  //  EVENTPREPARATION_HH
