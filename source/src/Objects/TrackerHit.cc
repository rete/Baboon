  /// \file TrackerHit.cc
/*
 *
 * TrackerHit.cc source template generated by fclass
 * Creation date : mer. juin 12 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Objects/TrackerHit.hh"

namespace baboon {

	TrackerHit::TrackerHit()
		: dEdx(0)
		,time(0)
		,eDeposit(0)
		,eDepositError(0) {

	}

	TrackerHit::~TrackerHit() {;}


	Return TrackerHit::SetPosition( const ThreeVector &pos ) {

		position = pos;
		return BABOON_SUCCESS();
	}

	Return TrackerHit::SetdEdx( const double &dedx ) {

		if( dedx < 0 )
			return BABOON_INVALID_PARAMETER("dE/dx < 0");
		dEdx = dedx;
		return BABOON_SUCCESS();
	}

	Return TrackerHit::SetTime( const double &t ) {

		if( t < 0 )
			return BABOON_INVALID_PARAMETER("Time < 0");
		time = t;
		return BABOON_SUCCESS();
	}

	Return TrackerHit::SetEDeposit( const double &e ) {

		if( e < 0 )
			return BABOON_INVALID_PARAMETER("Deposit energy < 0");
		eDeposit = e;
		return BABOON_SUCCESS();
	}

	Return TrackerHit::SetEDepositError( const double &eError ) {

		if( eError < 0 )
			return BABOON_INVALID_PARAMETER("Deposit energy error < 0");
		eDepositError = eError;
		return BABOON_SUCCESS();
	}

}  // namespace 
