  /// \file CoreToCoreAssociationAlgorithm.cc
/*
 *
 * CoreToCoreAssociationAlgorithm.cc source template generated by fclass
 * Creation date : ven. mai 24 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Algorithm/Association/CoreToCoreAssociationAlgorithm.hh"


namespace baboon {


	CoreToCoreAssociationAlgorithm::CoreToCoreAssociationAlgorithm()
		: AbstractAlgorithm("CoreToCoreAssociationAlgorithm") {
		needData = false;
	}


	CoreToCoreAssociationAlgorithm::~CoreToCoreAssociationAlgorithm() {

	}


	Return CoreToCoreAssociationAlgorithm::Init() {

		return BABOON_SUCCESS();
	}


	Return CoreToCoreAssociationAlgorithm::CheckConsistency() {

		return BABOON_SUCCESS();
	}

	Return CoreToCoreAssociationAlgorithm::Execute() {

		ShowerCollection *showerCollection = ShowerManager::GetInstance()->GetShowerCollection();

		if( showerCollection->size() <=1 )
			return BABOON_SUCCESS();

		/*

		for( unsigned int sh1=0 ; sh1<showerCollection->size() ; sh1++) {

			Shower *shower1 = showerCollection->at(sh1);
			CoreCollection *coreCollection1 = shower1->GetCoreCollection();

			for( unsigned int co1=0 ; co1<coreCollection1->size() ; co1++) {

				Core *core1 = coreCollection1->at(co1);

				for( unsigned int sh2=0 ; sh2<showerCollection->size() ; sh2++) {

					Shower *shower2 = showerCollection->at(sh2);

					if( shower1 == shower2 )
						continue;

					CoreCollection *coreCollection2 = shower2->GetCoreCollection();

					for( unsigned int co2=0 ; co2<coreCollection2->size() ; co2++) {

						Core *core2 = coreCollection2->at(co2);

						if( core1 == core2 )
							continue;

						if( this->CanAssociateCores( core1 , core2 ) ) ;

					}
				}
			}
		}
*/

		return BABOON_SUCCESS();
	}


	Return CoreToCoreAssociationAlgorithm::End() {

		return BABOON_SUCCESS();
	}

	bool CoreToCoreAssociationAlgorithm::CanAssociateCores( Core *core1 , Core *core2 ) {

		return false;
	}




}  // namespace 

