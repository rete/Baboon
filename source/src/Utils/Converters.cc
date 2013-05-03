/*
 *
 * Converters.cc cpp file template generated by fclass
 * Creation date : Wed Mar  6 11:09:36 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#include "Utils/Converters.hh"

using namespace std ;

namespace sdhcal {

	Converter::Converter() {}

	Converter::~Converter() {}

	vector<double> Converter::ThreeVectorToStdVector(const ThreeVector& g4Vec) {

		vector<double> vec;
		vec.push_back(g4Vec.x());
		vec.push_back(g4Vec.y());
		vec.push_back(g4Vec.z());

		return vec;
	}

	std::vector<double> Converter::StdStringVectorToStdDoubleVector(const std::vector<std::string>& v) {
		vector<double> vec;
		for(unsigned int i=0; i<v.size() ; i++) vec.push_back( atof( v.at(i).c_str() ) );
	}

	std::map<double,double>
	Converter::StdVectorsToStdMap(const std::vector<double>& v1 ,const std::vector<double>& v2 ) {
		if(v1.size() != v2.size()) throw domain_error("Tried to merge to vector with different sizes!");
		map<double,double> m;
		for(unsigned int i=0 ; i<v1.size() ; i++) m[v1.at(i)] = v2.at(i);
	}


	IMPL::CalorimeterHitImpl *
		Converter::CopyCalorimeterHitImpl(IMPL::CalorimeterHitImpl *hit) {
		IMPL::CalorimeterHitImpl *newHit = new IMPL::CalorimeterHitImpl();
		newHit->setCellID0( hit->getCellID0() );
		newHit->setCellID1( hit->getCellID1() );
		newHit->setEnergy( hit->getEnergy() );
		newHit->setEnergyError( hit->getEnergyError() );
		newHit->setTime( hit->getTime() );
		float pos[3]; for(int i=0;i<3;i++) pos[i] = hit->getPosition()[i];
		newHit->setPosition( pos );
		newHit->setType( hit->getType() );

		return newHit;
	}



}