/*
 *
 * AlgorithmManager.hh header template generated by fclass
 * Creation date : Thu Mar 14 22:21:50 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#ifndef ALGORITHMMANAGER_HH
#define ALGORITHMMANAGER_HH




#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <map>

// cfgparser includes
#include "CfgParser/Data.hh"
#include "CfgParser/CfgParser.hh"

// baboon includes
#include "Algorithm/AbstractAlgorithm.hh"
#include "Exception.hh"
#include "Utilities/Internal.hh"
#include "Utilities/ReturnValues.hh"


namespace baboon {


	typedef std::map< std::string , AbstractAlgorithm* > AlgorithmMap;

	/*!
	 * Class AlgorithmManager.
	 */

	class AlgorithmManager {

		protected :

			/*!
			 *
			 * Default Constructor
			 *
			 */
			AlgorithmManager();

			/*!
			 *
			 * Default Destructor. Not virtual because no inheritance
			 *
			 */
			~AlgorithmManager();

			static AlgorithmManager* instance;
			AlgorithmMap algorithmMap;
			std::string cfgFileName;

			/*!
			 *
			 * Load the algorithm parameters from the cfg file
			 *
			 */
			Return LoadAlgorithms();

			/*!
			 *
			 * @brief Print the algorithm header while initializing.
			 *
			 */
			void PrintAlgorithmHeader();

		public :

			/*!
			 *
			 * Only way to get an instance of this class. Part of the singleton patter
			 *
			 */
			static AlgorithmManager* GetInstance();

			/*!
			 *
			 * Only way to delete the object. Part of the singleton pattern
			 *
			 */
			static void Kill();

			/*!
			 *
			 * Register an algorithm
			 *
			 */
			Return RegisterAlgorithm( AbstractAlgorithm* );

			/*!
			 *
			 * Return true if the algorithm is registered.
			 * Useful to avoid exception while getting or registering an algorithm
			 *
			 */
			bool AlgorithmIsRegistered(const std::string&);

			/*!
			 *
			 * Return an algorithm which have already been registered via RegisterAlgorithm()
			 *
			 */
			AbstractAlgorithm* GetAlgorithm(const std::string&);

			/*!
			 *
			 * Set the cfg file name where the algorithm parameters are stored
			 *
			 */
			inline void SetConfigFileName(const std::string& cfgfName)
				{ cfgFileName = cfgfName; }

			/*!
			 *
			 * Initialize the beast.
			 *
			 */
			Return Initialize();


	};

}

#endif  // ALGORITHMMANAGER_HH
