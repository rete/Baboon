/*
 *
 * AbstractAlgorithm.hh header template generated by fclass
 * Creation date : Wed Mar  6 14:40:29 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#ifndef ABSTRACTALGORITHM_HH
#define ABSTRACTALGORITHM_HH

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <stdexcept>
#include <algorithm>

// cfgparser includes
#include "CfgParser/Data.hh"

// xml includes
#include "Xml/tinyxml.h"

// sdhcal includes
#include "Utilities/ReturnValues.hh"
#include "Reconstruction/Tag.hh"
#include "Exception.hh"

namespace baboon {

	/*!
	 * Class AbstractAlgorithm. Abstract class for Algorithm in SDHCAL studies. \n
	 * Define an interface for algorithm objects. Users have to implement their
	 * own algorithm via the four methods Init(), CheckConsistency(), Execute() and End(). \n
	 * Have to be passed to the AlgorithmManager while instantiating :\n
	 *\n
	 * @code
	 *
	 * AlgorithmManager *algoMan = AlgorithmManager::GetInstance();
	 * algoMan->SetConfigFileName("algoConfig.cfg");
	 * algoMan->RegisterAlgorithm( new MyAlgorithm() );
	 * ...
	 * MyAlgorithm *myAlgo;
	 * if( algoMan->AlgorithmIsRegistered("MyAlgorithm") ) {
	 *     myAlgo = (MyAlgorithm *) algoMan->GetAlgorithm("MyAlgorithm");
	 * }
	 *
	 * @endcode
	 */

	class AbstractAlgorithm {


		public :

			/*!
			 *
			 * @brief Default Constructor with algorithm name
			 *
			 */
			AbstractAlgorithm(const std::string& n)
				: name(n) , needData(false) {}

			/*!
			 *
			 * @brief Default Destructor
			 *
			 */
			virtual ~AbstractAlgorithm() {;}


			/*!
			 *
			 * @brief Method to be used by the user to start the algorithm
			 *
			 */
			void Process() {

				this->Init();
				Return ret = this->CheckConsistency();
				if( ret.OK ) {
					this->Execute();
					this->End();
				}
				else {
					std::cerr << ret.message << std::endl;
					throw AlgorithmException( ret.message );
				}

			}


			/*!
			 *
			 * @brief Return the algorithm name
			 *
			 */
			inline const std::string& GetName()
				{ return name; }


			/*!
			 *
			 * @brief Set the data for the algorithm
			 *
			 */
			inline void SetData(const cfgparser::Data& d)
				{ data = d; }

			/*!
			 *
			 * @brief Return the data set of the algorithm
			 *
			 */
			inline const cfgparser::Data& GetData()
				{ return data; }

			/*!
			 *
			 * @brief Return true if the algorithm need to be parameterized
			 *
			 */
			inline bool NeedSomeData()
				{ return needData; }


			/*!
			 *
			 * @brief Return true if the algorithm need to parameterized
			 *
			 */
			inline bool NeedParameters()
				{ return needParams; }


			/*!
			 *
			 * @brief Add a new tag to be avoided while running
			 *
			 */
			inline void AddAvoidedTag( Tag fTag ) {

				if( std::find( avoidedTags.begin(),avoidedTags.end(),fTag) == avoidedTags.end() )
					avoidedTags.push_back(fTag);
			}

		protected :

			/*!
			 *
			 * @brief
			 *
			 * */
//			virtual void Load( const TiXmlHandle *const pXmlHandle ) = 0;


			/*!
			 *
			 * @brief Initialize the algorithm, i.e by initializing specific variables
			 *
			 */
			virtual void Init() = 0;


			/*!
			 *
			 * @brief Execute the algorithm
			 *
			 */
			virtual void Execute() = 0;


			/*!
			 *
			 * @brief Finalize the algorithm
			 *
			 */
			virtual void End() = 0;


			/*!
			 *
			 * @brief Allow to check if everything is well set in the algorithm before starting it
			 *
			 */
			virtual Return CheckConsistency() = 0;


			std::string name;
			cfgparser::Data data;
			bool needData;
			bool needParams;
			std::vector<Tag> avoidedTags;

	};

}  // end of namespace

#endif  // ABSTRACTALGORITHM_HH
