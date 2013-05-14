  /// \file ShowerSplitterProcessor.cc
/*
 *
 * ShowerSplitterProcessor.cc source template generated by fclass
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


#include "OverlayEstimatorProcessor.hh"


OverlayEstimatorProcessor aOverlayEstimatorProcessor;


OverlayEstimatorProcessor::OverlayEstimatorProcessor()
	: marlin::Processor("OverlayEstimatorProcessor") {

	  _description = "OverlayEstimatorProcessor to estimate the overlap of two showers SDHCAL";

	  // register steering parameters: name, description, class-variable, default value

//	  std::vector<std::string> hcalCollections;
//	  hcalCollections.push_back(std::string("HCALBarrel"));
//	  registerInputCollections( LCIO::CALORIMETERHIT,
//				   "HCALCollections",
//				   "HCAL Collection Names",
//				   _hcalCollections,
//				   hcalCollections);
//
//	  registerOutputCollection( LCIO::LCRELATION,
//				    "RelationOutputCollection" ,
//				    "CaloHit Relation Collection" ,
//				    _outputRelCollection ,
//				    std::string("RelationCalorimeterHit")) ;
//
//	  registerProcessorParameter( "Energy" ,
//				      "Pion energy",
//				      energy_,
//				      (int) 0 );
//
//
//	  std::vector<float> thresholdHcal;
//	  thresholdHcal.push_back(0.114);
//	  thresholdHcal.push_back(1.45);
//	  thresholdHcal.push_back(3.80);
//	  registerProcessorParameter("HCALThreshold" ,
//	  			       "Threshold for HCAL Hits in GeV" ,
//	  			       _thresholdHcal,
//	  			       thresholdHcal);
//
//	  registerProcessorParameter( "RootFileName" ,
//				      "Name of the ROOT file where tree is stored",
//				      treeFileName_,
//				      std::string("showers.root") );
//
//	  std::vector<int> hcalLayers;
//	  hcalLayers.push_back(48);
//
//	  registerProcessorParameter("HCALLayers" ,
//				     "Index of HCal Layers" ,
//				     _hcalLayers,
//				     hcalLayers);

}

OverlayEstimatorProcessor::~ShowerSplitterProcessor() {

}


void OverlayEstimatorProcessor::processRunHeader( LCRunHeader* run) {

}

void OverlayEstimatorProcessor::processEvent( LCEvent * evt ) {

}


void OverlayEstimatorProcessor::check() {

}

void OverlayEstimatorProcessor::end() {

}

