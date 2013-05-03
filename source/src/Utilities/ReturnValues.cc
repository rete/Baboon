  /// \file ReturnValues.cc
/*
 *
 * ReturnValues.cc source template generated by fclass
 * Creation date : dim. avr. 28 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Utilities/ReturnValues.hh"

using namespace std;

//class Return;

namespace sdhcal {



	Return S_OK( const string &message ) {

		Return ret;
		ret.OK = true;
		ret.message = message;
		return ret;

	}

	Return S_ERROR( const string &message ) {

		Return ret;
		ret.OK = false;
		ret.message = message;
		return ret;

	}




}  // namespace 
