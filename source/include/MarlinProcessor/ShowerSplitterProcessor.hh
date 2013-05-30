  /// \file ShowerSplitterProcessor.hh
/*
 *
 * ShowerSplitterProcessor.hh header template generated by fclass
 * Creation date : jeu. mai 9 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef SHOWERSPLITTERPROCESSOR_HH
#define SHOWERSPLITTERPROCESSOR_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 


#include "MarlinProcessor/BaboonProcessor.hh"
#include "Reconstruction/SimpleEnergyCalculator.hh"

/* 
 * Class ShowerSplitterProcessor
 * Inherits from base class BaboonProcessor
 */ 

class ShowerSplitterProcessor : public baboon::BaboonProcessor {

	public:

		/*! Default Constructor */
		ShowerSplitterProcessor();

		/*! Default Destructor */
		virtual ~ShowerSplitterProcessor();


		virtual Processor *newProcessor()
			{ return new ShowerSplitterProcessor(); }

		/*!
		 *
		 * @brief Must be defined by the user. Init the BABOON processor
		 *
		 */
		virtual baboon::Return Init();

		/*!
		 *
		 * @brief Must be defined by the user. Process the run header.
		 *
		 */
		virtual baboon::Return ProcessRunHeader( EVENT::LCRunHeader* run  );

		/*!
		 *
		 * @brief Must be defined by the user. Process an event in the BABOON framework
		 *
		 */
		virtual baboon::Return ProcessEvent( const unsigned int &evtNb );

		/*!
		 *
		 * @brief Must be defined by the user. Check the event.
		 *
		 */
		virtual baboon::Return Check( EVENT::LCEvent * evt );

		/*!
		 *
		 * @brief Must be defined by the user. Called after processing all the events.
		 *
		 */
		virtual baboon::Return End();


};  // class 


#endif  //  SHOWERSPLITTERPROCESSOR_HH
