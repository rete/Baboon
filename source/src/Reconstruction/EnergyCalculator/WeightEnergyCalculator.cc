  /// \file WeightEnergyCalculator.cc
/*
 *
 * WeightEnergyCalculator.cc source template generated by fclass
 * Creation date : mer. mai 29 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Reconstruction/EnergyCalculator/WeightEnergyCalculator.hh"

using namespace std;

namespace baboon {

	WeightEnergyCalculator::WeightEnergyCalculator()
		: EnergyCalculator() {

		params = new double[9];
		params[0] = 0.0227148;
		params[1] = 4.77839e-05;
		params[2] = -4.05221e-08;
		params[3] = 0.112287;
		params[4] = 6.24094e-06;
		params[5] = 2.9642e-08;
		params[6] = 0.156503;
		params[7] = 0.00010774;
		params[8] = 4.97735e-08;

		shower = 0;
	}

	WeightEnergyCalculator::~WeightEnergyCalculator() {

		delete [] params;
	}


	Return WeightEnergyCalculator::SetShower( Shower *sh ) {

		if( sh == 0 )
			return BABOON_INVALID_PARAMETER("Assertion shower != 0 failed");

		shower = sh;
		return BABOON_SUCCESS();
	}



	Return WeightEnergyCalculator::CalculateEnergy() {

		if( shower == 0 )
			return BABOON_NOT_INITIALIZED("Can't compute energy without setting a shower before.");

		const DoubleVector& hitWeights = shower->GetHitWeights();
		CaloHitCollection *caloHitCollection = shower->GetCaloHitCollection();

		energy = 0;

		if( caloHitCollection->empty() )
			return BABOON_SUCCESS();

		double nbOfHits = 0;
		double nbThreshold1 = 0;
		double nbThreshold2 = 0;
		double nbThreshold3 = 0;

		for( unsigned int h=0 ; h<caloHitCollection->size() ; h++ ) {

			nbOfHits += hitWeights.at(h);
			CaloHit* caloHit = caloHitCollection->at(h);
			if( caloHit->GetThreshold() == fCaloHitThr1 ) nbThreshold1 += hitWeights.at(h);
			else if( caloHit->GetThreshold() == fCaloHitThr2 ) nbThreshold2 += hitWeights.at(h);
			else if( caloHit->GetThreshold() == fCaloHitThr3 ) nbThreshold3 += hitWeights.at(h);

		}

		double alpha = params[0] + params[1]*( nbOfHits ) + params[2]*( nbOfHits*nbOfHits );
		double beta = params[3] + params[4]*( nbOfHits ) + params[5]*( nbOfHits*nbOfHits );
		double gamma = params[6] + params[7]*( nbOfHits ) + params[8]*( nbOfHits*nbOfHits );

		energy = alpha*nbThreshold1 + beta*nbThreshold2 + gamma*nbThreshold3;
		return BABOON_SUCCESS();
	}

	Return WeightEnergyCalculator::SetParameter( const int &i , const double &value ) {

		if( i > 8 )
			return BABOON_OUT_OF_RANGE("Out of range for params[i]");
		else params[i] = value;
		return BABOON_SUCCESS();
	}


}  // namespace 

