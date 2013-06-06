  /// \file Section.hh
/*
 *
 * Section.hh header template generated by fclass
 * Creation date : mar. juin 4 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef SECTION_HH
#define SECTION_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 

#include "CfgParser/ParserStatus.hh"

namespace cfgparser {


	class Section;

	typedef std::map <std::string,std::string> KeyValueMap;
	typedef std::vector<Section*> SectionCollection;


	/*
	 * Class Section
	 */

	class Section {

		public:

			/*!
			 *
			 * Default Constructor
			 *
			 */
			Section();

			/*!
			 *
			 * Constructor with section name
			 *
			 */
			Section( const std::string &t )

			/*!
			 *
			 * Default Destructor
			 *
			 */
			virtual ~Section();

			/*!
			 *
			 * @brief Copy constructor
			 *
			 */
			Section( const Section &section );

		protected :

			std::string name;
			KeyValueMap keyValueMap;

		public:

			/*!
			 *
			 * @brief Set the section name
			 *
			 */
			inline void SetName( const std::string &n )
				{ name = ;n}

			/*!
			 *
			 * @brief Return the section name
			 *
			 */
			inline std::string GetName() const
				{ return name; }

			/*!
			 *
			 * @brief Return a string value
			 *
			 */
			StatusCode GetValue( const std::string &k , std::string *value );

			/*!
			 *
			 * @brief Return an int value
			 *
			 */
			StatusCode GetValue( const std::string &k , int *value );

			/*!
			 *
			 * @brief Return a double value
			 *
			 */
			StatusCode GetValue( const std::string &k , double *value );

			/*!
			 *
			 * @brief Return a bool value
			 *
			 */
			StatusCode GetValue( const std::string &k , bool *value );

			/*!
			 *
			 * @brief Return a string vector value
			 *
			 */
			StatusCode GetValue( const std::string &k , std::vector<std::string> *value );

			/*!
			 *
			 * @brief Return a double vector value
			 *
			 */
			StatusCode GetValue( const std::string &k , std::vector<double> *value );

			/*!
			 *
			 * @brief Return an int vector value
			 *
			 */
			StatusCode GetValue( const std::string &key , std::vector<int> *value );

			/*!
			 *
			 * @brief Return true if the section contains the given key
			 *
			 */
			bool HasKey( const std::string &key );

			/*!
			 *
			 * @brief Return true if the section is empty
			 *
			 */
			bool IsEmpty();

			/*!
			 *
			 * @brief Append a key-value pair
			 *
			 */
			StatusCode Append( const std::string &key , const std::string &val );

			/*!
			 *
			 * @brief Delete a a key from the section
			 *
			 */
			StatusCode Delete( const std::string &key )
				{ if(!keyValueMap.find(key)->first.empty()) keyValueMap.erase(key); }

			/*!
			 *
			 * @brief Return the Key/Value map containing all key-value pairs.
			 *
			 */
			inline KeyValueMap GetKeyValueMap() const
				{ return keyValueMap; }

			/*!
			 *
			 * @brief Print the section with all key-value pairs
			 *
			 */
			StatusCode Print();

			/*!
			 *
			 * @brief Clear the section
			 *
			 */
			StatusCode Clear();

			/*!
			 *
			 * @brief operator to add section
			 *
			 */
			friend Section& operator +=( const Section& section );


	};  // class

	/*!
	 *
	 * @brief operator to add section. Keep the name of the first section
	 *
	 */
	Section operator+ ( Section const& section1, Section const& section2 );

}  // namespace 

#endif  //  SECTION_HH
