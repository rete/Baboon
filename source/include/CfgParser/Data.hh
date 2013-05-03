/*
 *
 * //Data storage for configuration file parser
 *
 * Use and copying of this software and preparation of derivative works
 * based upon this software are permitted. Any copy of this software or
 * of any derivative work must include this paragraph.
 *
 * Written by : rete
 *
 */

#ifndef DATA_HH
#define DATA_HH

#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>
#include <list>
#include <typeinfo>


// cfgparser includes
#include "CfgParser/CfgParserException.hh"
#include "CfgParser/StrOp.hh"


namespace cfgparser {


	class Data;

	typedef std::map <std::string,std::string> KeyValueMap;
	typedef std::vector<Data> DataList;

	/*!
	 *
	 * Data class.
	 * Data is identified by its type.
	 * Contains a map of key and value.
	 *
	 */
	class Data {


		protected :

			std::string type;
			KeyValueMap keyValueMap;

		public :


			/*!
			 *
			 * @brief Default constructor
			 *
			 */
			Data();

			/*!
			 *
			 * @brief Explicit contructor with data type name
			 *
			 */
			explicit Data( const std::string &t );

			/*!
			 *
			 * @brief Copy constructor
			 *
			 */
			Data( const Data &data );

			/*!
			 *
			 * @brief Default destructor
			 *
			 */
			~Data();

			/*!
			 *
			 * @brief Set the data type
			 *
			 */
			inline void SetType( const std::string &t )
				{ type = t;}

			/*!
			 *
			 * @brief Return the data type
			 *
			 */
			inline std::string GetType() const
				{ return type; }

			/*!
			 *
			 * @brief Return a string value with key 'k'
			 *
			 */
			void GetValue( const std::string &k , std::string *v);

			/*!
			 *
			 * @brief Return an int value with key 'k'
			 *
			 */
			void GetValue( const std::string &k , int *v);

			/*!
			 *
			 * @brief Return a double value with key 'k'
			 *
			 */
			void GetValue( const std::string &k , double *v);

			/*!
			 *
			 * @brief Return a bool value with key 'k'
			 *
			 */
			void GetValue( const std::string &k , bool *v);

			/*!
			 *
			 * @brief Return a string vector value with key 'k'
			 *
			 */
			void GetValue( const std::string &k , std::vector<std::string> *v);

			/*!
			 *
			 * @brief Return a double vector value with key 'k'
			 *
			 */
			void GetValue( const std::string &k , std::vector<double> *v );

			/*!
			 *
			 * @brief Return an int vector value with key 'k'
			 *
			 */
			void GetValue( const std::string &k , std::vector<int> *v );

			/*!
			 *
			 * @brief Return the Key/Value map containing all key-value pairs.
			 *
			 */
			inline KeyValueMap GetKeyValueMap() const
				{ return keyValueMap; }

			/*!
			 *
			 * @brief Return true if the data contains a key named 'k'
			 *
			 */
			bool HasKey( const std::string &k );

			/*!
			 *
			 * @brief Return true if the data is empty
			 *
			 */
			bool IsEmpty();

			/*!
			 *
			 * @brief Append a key-value pair with repective name 'k' and 'v'
			 *
			 */
			void Append( const std::string &k , const std::string &v );

			/*!
			 *
			 * @brief Delete a a key from the data block
			 *
			 */
			void DeleteKey( const std::string &k )
				{ if(!keyValueMap.find(k)->first.empty()) keyValueMap.erase(k); }

			/*!
			 *
			 * @brief Print the data with all key-value pair
			 *
			 */
			void PrintData();

			/*!
			 *
			 * @brief Clear the data block
			 *
			 */
			void Clear();

			/*!
			 *
			 * @brief operator to add data
			 *
			 */
			Data& operator +=(const Data& data);

	};
	
	/*!
	 *
	 * @brief operator to add data
	 *
	 */
	Data operator+ (Data const& data1, Data const& data2);




}  //  end namespace cfgparser


#endif  // DATA_HH
