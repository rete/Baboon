/*
 *
 * AlgorithmManager.cc cpp file template generated by fclass
 * Creation date : Thu Mar 14 22:21:50 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#include "Managers/AlgorithmManager.hh"



using namespace std;
using namespace cfgparser;



namespace baboon {


	AlgorithmManager* AlgorithmManager::instance = NULL;



	AlgorithmManager::AlgorithmManager() {}



	AlgorithmManager::~AlgorithmManager() {

		AlgorithmMap::iterator it;
		for( it=algorithmMap.begin() ; it!=algorithmMap.end() ; it++ ) {
			if( (*it).second != NULL ) delete (*it).second;
		}
		algorithmMap.clear();
	}


	AlgorithmManager* AlgorithmManager::GetInstance() {
		if(instance == NULL) instance = new AlgorithmManager();
		return instance;
	}

	void AlgorithmManager::Kill() {
		if(instance != NULL) {
			delete instance;
			instance = NULL;
		}
	}

	bool AlgorithmManager::AlgorithmIsRegistered(const string& algoName) {

		AlgorithmMap::iterator it;
		for( it=algorithmMap.begin() ; it!=algorithmMap.end() ; it++ ) {
			if( (*it).first == algoName ) return true;
		}
		return false;
	}

	void AlgorithmManager::RegisterAlgorithm( AbstractAlgorithm* algo)
		throw ( AlgorithmException , Exception ) {

		string algoName = algo->GetName();

		if( !this->AlgorithmIsRegistered( algoName ) ) {
			algorithmMap[algoName] = algo;
		}
		else {
			delete algo;
			throw AlgorithmException("Algorithm "+algoName+" is already registered!");
		}
	}

	AbstractAlgorithm* AlgorithmManager::GetAlgorithm( const std::string& algoName )
		throw ( AlgorithmException , Exception ) {

		AlgorithmMap::iterator it;
		AbstractAlgorithm *algo;
		bool found = false;

		for( it=algorithmMap.begin() ; it!=algorithmMap.end() ; it++ ) {
			if( (*it).first == algoName ) {
				found = true;
				algo = (*it).second;
			}
		}
		if( !found ) {
			throw AlgorithmException(
					"Algorithm "+algoName+" doesn't exists or is not registered in AlgorithmManager" );
		}
		else return algo;
	}

	void AlgorithmManager::LoadAlgorithms() {

		AlgorithmMap::iterator it;

		CfgParser *parser = new CfgParser;
		parser->SetConfigFileName(cfgFileName);
		parser->Read();

		for( it=algorithmMap.begin() ; it!=algorithmMap.end() ; it++ ) {
			string algoName = (*it).first;

			try {

				Data algoData = parser->GetData(algoName);
				(*it).second->SetData(algoData);

			} catch ( std::exception &e ) {

				if( it->second->NeedSomeData() ) {
					throw AlgorithmException( "Algorithm " + algoName + " need some data to be processed. "
							"Please add a section in " + cfgFileName);
				}
				else continue;
			}
		}

		delete parser;
	}


	void AlgorithmManager::PrintAlgorithmHeader() {

		AlgorithmMap::iterator it;

		cout << "*********************************************************" << endl;
		cout << "* Algorithm Manager registered the following algorithms :" << endl;
		cout << "*" << endl;
		for( it=algorithmMap.begin() ; it!=algorithmMap.end() ; it++ ) {
			string algoName = (*it).first;
			cout << "*    - " << algoName << endl;
		}
		cout << "*" << endl;
		cout << "*********************************************************" << endl;

	}


	void AlgorithmManager::Initialize() {

		instance->LoadAlgorithms();
		instance->PrintAlgorithmHeader();
	}


}
