/// \file CfgParser.hh
 /*
 *
 * //Generic configuration file parser
 *
 * Use and copying of this software and preparation of derivative works
 * based upon this software are permitted. Any copy of this software or
 * of any derivative work must include this paragraph.
 *
 * @author rete
 *
 */

#ifndef CFGPARSER_HH
#define CFGPARSER_HH

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <map>
#include <vector>


// cfgparser includes
#include "CfgParser/Data.hh"
#include "CfgParser/CfgParserException.hh"
#include "CfgParser/StrOp.hh"


namespace cfgparser {


	/*!
	 *
	 * Configuration file Parser.
	 * Read a configuration file in a specific syntax.
	 * Store all data given in Data objects.
	 * See getter functions to know how to use it.
	 *
	 */

	class CfgParser {

		public :

			/*!
			 *
			 * @brief Default constructor
			 *
			 */
			CfgParser();

			/*!
			 *
			 * @brief Constructor with file name
			 *
			 */
			CfgParser( const std::string &fileName );

			/*!
			 *
			 * @brief Default destructor
			 *
			 */
			~CfgParser();

			/*!
			 *
			 * @brief read the cfg file. Throw if the file name is not set or can't be read.
			 *
			 */
			void Read() throw ( CfgParserException , std::exception );

			/*!
			 *
			 * @brief return a data block in the cfg file
			 *
			 */
			Data GetData( const std::string &t );

			/*!
			 *
			 * @brief return the cfg filename
			 *
			 */
			const std::string &GetFileName()
				{ return configFileName; }

			/*!
			 *
			 * @brief set the cfg file name
			 *
			 */
			inline void SetConfigFileName(const std::string& fileName)
				{ configFileName = fileName; }

			/*!
			 *
			 * @brief return a string value in the data type 't' with key 'k'.
			 *
			 */
			void GetValue( const std::string &t , const std::string &k , std::string *v);

			/*!
			 *
			 * @brief return an int value in the data type 't' with key 'k'.
			 *
			 */
			void GetValue( const std::string &t , const std::string &k , int *v);

			/*!
			 *
			 * @brief return a double value in the data type 't' with key 'k'.
			 *
			 */
			void GetValue( const std::string &t , const std::string &k , double *v);

			/*!
			 *
			 * @brief return a bool value in the data type 't' with key 'k'.
			 *
			 */
			void GetValue( const std::string &t , const std::string &k , bool *v);

			/*!
			 *
			 * @brief return a string vector value in the data type 't' with key 'k'.
			 *
			 */
			void GetValue( const std::string &t , const std::string &k , std::vector<std::string> *v);

			/*!
			 *
			 * @brief return a double vector value in the data type 't' with key 'k'.
			 *
			 */
			void GetValue( const std::string &t , const std::string &k , std::vector<double> *v);

			/*!
			 *
			 * @brief return an int vector value in the data type 't' with key 'k'.
			 *
			 */
			void GetValue( const std::string &t , const std::string &k , std::vector<int> *v);

			/*!
			 *
			 * @brief return true if the file is already read.
			 *
			 */
			bool IsAlreadyRead()
				{ return alreadyRead; }

			/*!
			 *
			 * @brief print the data block list (should be used after reading).
			 *
			 */
			void PrintDataList();

			/*!
			 *
			 * @brief return the total data list content in the cfg file (should be used after reading).
			 *
			 */
			DataList GetDataList()
				{ return dataList; }


		protected :

			std::fstream *file;
			std::string configFileName;
			DataList dataList;
			std::map <int,std::string> idLine;
			std::string optionRead;
			bool alreadyRead;

	};

}  //  end namespace cfgparser


#endif  // CFGPARSER_HH
