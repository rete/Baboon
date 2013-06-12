  /// \file EventDisplay.hh
/*
 *
 * EventDisplay.hh header template generated by fclass
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


#ifndef EVENTDISPLAY_HH
#define EVENTDISPLAY_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 


#include "MarlinProcessor/BaboonProcessor.hh"

/* 
 * Class EventDisplay
 * Inherits from base class BaboonProcessor
 */ 

class EventDisplay : public baboon::BaboonProcessor {

	public:

		/*! Default Constructor */
		EventDisplay();

		/*! Default Destructor */
		virtual ~EventDisplay();


		virtual Processor *newProcessor()
			{ return new EventDisplay(); }

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
		virtual baboon::Return ProcessEvent( EVENT::LCEvent * );

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

	protected:

		std::string displayMode;


};  // class 


#endif  //  EVENTDISPLAY_HH
