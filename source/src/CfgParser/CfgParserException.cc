/*
 *
 * //Parser exception classes for CfgParser
 *
 * Use and copying of this software and preparation of derivative works
 * based upon this software are permitted. Any copy of this software or
 * of any derivative work must include this paragraph.
 *
 * Written by : rete
 *
 */

using namespace std;



#include "CfgParser/CfgParserException.hh"



namespace cfgparser {


	const char* CfgParserException::what() const throw() {return error.c_str();}

	const char* DataException::what() const throw() {return error.c_str();}


}  //  end namespace cfgparser
