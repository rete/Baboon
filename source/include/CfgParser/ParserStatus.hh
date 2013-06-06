  /// \file ParserStatus.hh
/*
 *
 * ParserStatus.hh header template generated by fclass
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


#ifndef PARSERSTATUS_HH
#define PARSERSTATUS_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 

#include "CfgParser/CfgParserException.hh"

#define CFGPARSER_THROW_RESULT_IF( Status , Operator , Command )                                    \
	{                        																		 \
                                                     	 	 	 	 	 	 	 	 	 	 	 	 \
		if( Status.fParserStatus Operator Command.fParserStatus ) { 								 \
																									 \
	        std::cerr << #Command << " throw " << Command.ToString() << std::endl;                   \
	        std::cerr << "    in function : " << __FUNCTION__ << std::endl							 \
	        std::cerr << "    in file :     " << __FILE__ << " line#: " << __LINE__ << std::endl;    \
	        std::cerr << "    message :     " << Command.message << std::endl;                       \
	        throw CfgParserException(Command.message);												 \
		}																							 \
	}

namespace cfgparser {


	enum ParserStatus {
		fParserSuccess,
		fParserNoSectionError,
		fParserDuplicateSectionError,
		fParserNoOptionError,
		fParserInterpolationError,
		fParserInterploationDepthError,
		fParserInterpolationMissingOptionError,
		fParserInterpolationSyntaxError,
		fParserMissingSectionNameError,
		fParserInvalidSectionKey,
		fParserParsingError
	};



	class StatusCode {

		ParserStatus fParserStatus;
		std::string message;
		std::string ToString();
	};

	StatusCode CFGPARSER_SUCCESS();

	StatusCode CFGPARSER_NO_SECTION_ERROR();

	StatusCode CFGPARSER_DUPLICATE_SECTION();

	StatusCode CFGPARSER_NO_OPTION_ERROR();

	StatusCode CFGPARSER_INTERPOLATION_ERROR();

	StatusCode CFGPARSER_INTERPOLATION_DEPTH_ERROR();

	StatusCode CFGPARSER_INTERPOLATION_MISSING_OPTION_ERROR();

	StatusCode CFGPARSER_INTERPOLATION_SYNTAX_ERROR();

	StatusCode CFGPARSER_SECTION_HEADER_ERROR();

	StatusCode CFGPARSER_INVALID_SECTION_KEY();

	StatusCode CFGPARSER_MISSING_SECTION_NAME_ERROR();

	StatusCode CFGPARSER_PARSING_ERROR();

}  // namespace 

#endif  //  PARSERSTATUS_HH
