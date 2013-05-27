  /// \file GShowerSplitterProcessor.hh
/*
 *
 * GShowerSplitterProcessor.hh header template generated by fclass
 * Creation date : lun. mai 20 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef GSHOWERSPLITTERPROCESSOR_HH
#define GSHOWERSPLITTERPROCESSOR_HH

#include "GraphicalProcessor/GraphicalProcessor.hh"

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 

/* 
 * Class GShowerSplitterProcessor
 * Inherits from base class GraphicalProcessor
 */ 

class GShowerSplitterProcessor : public baboon::GraphicalProcessor {

	public:
		/*! Default Constructor */
		GShowerSplitterProcessor();
		/*! Default Destructor */
		virtual ~GShowerSplitterProcessor();

		Processor *newProcessor()
			{ return new GShowerSplitterProcessor(); }

	protected:

		virtual baboon::Return Init();

		virtual baboon::Return ProcessRunHeader( EVENT::LCRunHeader* run  );

		virtual baboon::Return ProcessEvent();

		virtual baboon::Return Check( EVENT::LCEvent * evt );

		virtual baboon::Return End();


};  // class 


#endif  //  GSHOWERSPLITTERPROCESSOR_HH