  /// \file SdhcalConfig.cc
/*
 *
 * SdhcalConfig.cc source template generated by fclass
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



#include "Config/SdhcalConfig.hh"

using namespace std;
using namespace cfgparser;

namespace baboon {



	SdhcalConfig *SdhcalConfig::instance = 0;



	SdhcalConfig::SdhcalConfig() {
		parser = new CfgParser();
	}



	SdhcalConfig::~SdhcalConfig() {
		delete parser;
	}



	SdhcalConfig *SdhcalConfig::GetInstance() {

		if( instance == 0 ) instance = new SdhcalConfig();
		return instance;

	}



	void SdhcalConfig::Kill() {

		if( instance != 0 ) {
			delete instance;
			instance = 0;
		}
	}



	void SdhcalConfig::LoadFile( const std::string &cfgfName ) {

		cfgFileName = cfgfName;
		parser->SetConfigFileName(cfgFileName);
		parser->Read();

	}


	Data SdhcalConfig::GetData( const string &data ) {
		return parser->GetData(data);
	}

}  // namespace 

