  /// \file Tag.cc
/*
 *
 * Tag.cc source template generated by fclass
 * Creation date : lun. avr. 22 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * @author : rete
 */

#include "Reconstruction/Tag.hh"


namespace baboon {



	Tag StringToTag( const std::string tagString ) {

		if( tagString == "fIsolated" )             return fIsolated;
		else if( tagString == "fTrack" )           return fTrack;
		else if( tagString == "fTrackExtremity" )  return fTrackExtremity;
		else if( tagString == "fCore" )            return fCore;
		else if( tagString == "fNoise" )           return fNoise;
		else                                       return fUndefined;
	}

	std::string TagToString( Tag tag ) {

		      if( tag == fIsolated )       return std::string("fIsolated");
		else if( tag == fTrack )          return std::string("fTrack");
		else if( tag == fTrackExtremity ) return std::string("fTrackExtremity");
		else if( tag == fCore )           return std::string("fCore");
		else if( tag == fNoise )       return std::string("fNoise");
		else                               return std::string("fUndefined");
	}


//	BaseTag::BaseTag( const std::string &name )


}  // namespace

