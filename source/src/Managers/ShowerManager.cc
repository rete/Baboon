  /// \file ShowerManager.cc
/*
 *
 * ShowerManager.cc source template generated by fclass
 * Creation date : mar. mai 14 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Managers/ShowerManager.hh"

using namespace std;

namespace baboon {

	ShowerManager *ShowerManager::instance = 0;


	ShowerManager::ShowerManager() {}


	ShowerManager::~ShowerManager() {}


	ShowerManager *ShowerManager::GetInstance() {

		if( instance == 0 )
			instance = new ShowerManager();
		return instance;
	}

	void ShowerManager::Kill() {

		if( instance != 0 ) {
			delete instance;
			instance = 0;
		}
	}


	Return ShowerManager::AddShower( Shower *shower ) {

		if( showerCollection->empty() ) {
			showerCollection->push_back( shower );
			return S_OK();
		}
		ShowerCollection::iterator showerIt = std::find( showerCollection->begin() , showerCollection->end() , shower );

		if( showerIt == showerCollection->end() ) {
			showerCollection->push_back( shower );
			return S_OK();
		}

		return S_ERROR("While adding shower. Shower already registered by the Shower Manager!");
	}



	Return ShowerManager::RemoveShower( Shower *shower ) {

		if( showerCollection->empty() )
			return S_ERROR("While removing a shower. Shower collection is empty!");

		ShowerCollection::iterator showerIt = std::find( showerCollection->begin() , showerCollection->end() , shower );

		if( showerIt != showerCollection->end() ) {
			delete shower;
			showerCollection->erase( showerIt );
			return S_OK("Shower correctly removed");
		}

		return S_ERROR("While removing a shower. Shower was not registered by the Shower Manager!");
	}


	bool ShowerManager::ShowerContainsHit( Shower *shower , Hit *hit ) {

		return shower->Contains( hit );
	}


}  // namespace 
