  /// \file EventDisplay.cc
/*
 *
 * EventDisplay.cc source template generated by fclass
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


#include "MarlinProcessor/EventDisplay.hh"


EventDisplay aEventDisplay;


using namespace std;
using namespace baboon;
using namespace EVENT;

EventDisplay::EventDisplay()
	: BaboonProcessor("EventDisplay") {

	  _description = "EventDisplay to display event in SDHCAL";

	  registerProcessorParameter("displayMode" ,
			  	  	  "display mode : threshold (default), overlay" ,
			  	  	  displayMode,
			  	  	  string("threshold"));

}

EventDisplay::~EventDisplay() {}


Return EventDisplay::Init() {

	return BABOON_SUCCESS();
}

Return EventDisplay::ProcessRunHeader( LCRunHeader* run ) {

	return BABOON_SUCCESS();
}

Return EventDisplay::ProcessEvent( EVENT::LCEvent * evt ) {

	if( displayMode != "threshold" && displayMode != "overlay" )
		displayMode = "threshold";

	HitCollection *hitCollection = hitManager->GetHitCollection();

	if( hitCollection->empty() )
		return BABOON_NOT_INITIALIZED("No hit in the current event. Framework not initialized!");

	for(unsigned int j=0 ; j<hitCollection->size() ; j++) {

		IntVector ijk = hitCollection->at(j)->GetIJK();

		if( displayMode == "threshold" ) {

			if( hitCollection->at(j)->GetThreshold() == fThreshold1 ) {
				if( graphicalEnvironment )
					calorimeter->GetNodeAt(ijk.at(0),ijk.at(1),ijk.at(2))->GetVolume()->SetLineColor(kGreen);
			}
			if( hitCollection->at(j)->GetThreshold() == fThreshold2 ) {
				if( graphicalEnvironment )
					calorimeter->GetNodeAt(ijk.at(0),ijk.at(1),ijk.at(2))->GetVolume()->SetLineColor(kBlue);
			}
			if( hitCollection->at(j)->GetThreshold() == fThreshold3 ) {
				if( graphicalEnvironment )
					calorimeter->GetNodeAt(ijk.at(0),ijk.at(1),ijk.at(2))->GetVolume()->SetLineColor(kRed);
			}
		}
		else if( displayMode == "overlay" ) {
			if( hitCollection->at(j)->GetType() == 1 ) {
				if( graphicalEnvironment )
					calorimeter->GetNodeAt(ijk.at(0),ijk.at(1),ijk.at(2))->GetVolume()->SetLineColor(kGreen);
			}
			else if( hitCollection->at(j)->GetType() == 2 ) {
				if( graphicalEnvironment )
					calorimeter->GetNodeAt(ijk.at(0),ijk.at(1),ijk.at(2))->GetVolume()->SetLineColor(kBlue);
			}
			else if( hitCollection->at(j)->GetType() == 3 ) {
				if( graphicalEnvironment )
					calorimeter->GetNodeAt(ijk.at(0),ijk.at(1),ijk.at(2))->GetVolume()->SetLineColor(kRed);
			}
		}
	}

	return BABOON_SUCCESS();

}


Return EventDisplay::Check( LCEvent *evt ) {

	return BABOON_SUCCESS();
}

Return EventDisplay::End() {

	return BABOON_SUCCESS();
}


