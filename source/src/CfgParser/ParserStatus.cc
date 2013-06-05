  /// \file ParserStatus.cc
/*
 *
 * ParserStatus.cc source template generated by fclass
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


#include "ParserStatus.hh"

namespace cfgparser {


StatusCode ParserStatus::ToString() {

	if( fStatus == fParserSuccess )    return "fParserSuccess";
	else if( fStatus == fParserNoSectionError ) return "fParserNoSectionError";
	else if( fStatus == fParserDuplicateSectionError ) return "fParserDuplicateSectionError";
	else if( fStatus == fParserNoOptionError ) return "fParserNoOptionError";
	else if( fStatus == fParserInterpolationError ) return "fParserInterpolationError";
	else if( fStatus == fParserInterploationDepthError ) return "fParserInterploationDepthError";
	else if( fStatus == fParserInterploationDepthError ) return "fParserInterploationDepthError";
	else if( fStatus == fParserInterpolationMissingOptionError ) return "fParserInterpolationMissingOptionError";
	else if( fStatus == fParserInterpolationSyntaxError ) return "fParserInterpolationSyntaxError";
	else if( fStatus == fParserMissingSectionHeaderError ) return "fParserMissingSectionHeaderError";
	else if( fStatus == fParserInvalidSectionKey ) return "fParserInvalidSectionKey";
	else return "fParserParsingError";
}


StatusCode CFGPARSER_SUCCESS( std::string &message ) {
	StatusCode statusCode;
	statusCode.fStatus = fParserSuccess;
	statusCode.message = message;
	return statusCode;
}

StatusCode CFGPARSER_INVALID_SECTION_KEY( std::string &message ) {
	StatusCode statusCode;
	statusCode.fStatus = fParserInvalidSectionKey;
	statusCode.message = message;
	return statusCode;
}

StatusCode CFGPARSER_NO_SECTION_ERROR( std::string &message );

StatusCode CFGPARSER_DUPLICATE_SECTION( std::string &message );

StatusCode CFGPARSER_NO_OPTION_ERROR( std::string &message );

StatusCode CFGPARSER_INTERPOLATION_ERROR( std::string &message );

StatusCode CFGPARSER_INTERPOLATION_DEPTH_ERROR( std::string &message );

StatusCode CFGPARSER_INTERPOLATION_MISSING_OPTION_ERROR( std::string &message );

StatusCode CFGPARSER_INTERPOLATION_SYNTAX_ERROR( std::string &message );

StatusCode CFGPARSER_SECTION_HEADER_ERROR( std::string &message );

StatusCode CFGPARSER_PARSING_ERROR( std::string &message );


}  // namespace 

