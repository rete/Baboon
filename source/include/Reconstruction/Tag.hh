  /// \file Tag.hh
/*
 *
 * Tag.hh header template generated by fclass
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


#ifndef TAG_HH
#define TAG_HH


#include <string>


namespace sdhcal {

	/*!
	 * Enumerator for hit tagging. Allow to classify hits in SDHCAL.
	 */
	enum Tag {
		fIsolated,
		fTrackSegment,
		fTrack,
		fTrackExtremity,
		fCore,
		fCoreEdge,
		fUndefined
	};


	Tag StringToTag( const std::string tagString );

	std::string TagToString( Tag tag );


}  // namespace 

#endif  //  TAG_HH
