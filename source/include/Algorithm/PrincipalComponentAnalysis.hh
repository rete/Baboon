  /// \file PrincipalComponentAnalysis.hh
/*
 *
 * PrincipalComponentAnalysis.hh header template generated by fclass
 * Creation date : dim. avr. 21 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef PRINCIPALCOMPONENTANALYSIS_HH
#define PRINCIPALCOMPONENTANALYSIS_HH

#include "Algorithm/AbstractAlgorithm.hh"

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <map>
#include <numeric>

#include "TMatrixT.h"
#include "TMatrixDEigen.h"
#include "TVectorD.h"
#include "TPrincipal.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"

#include "Geometry/ThreeVector.hh"
#include "Objects/Hit.hh"
#include "Exception.hh"


namespace baboon {

	/*
	 * Class PrincipalComponentAnalysis
	 * Inherits from base class AbstractAlgorithm
	 */

	class PrincipalComponentAnalysis : public AbstractAlgorithm {

		public:

			/*! Default Constructor */
			PrincipalComponentAnalysis();

			/*! Default Destructor */
			virtual ~PrincipalComponentAnalysis();

			/**/
			inline void SetHitCollection( HitCollection *hitCol )
				{ hitCollection = hitCol; }

			inline TMatrixD GetEigenVectors()
				{ return eigenVectors; }

			inline TVectorD GetEigenValues()
				{ return eigenValues; }


		protected:

			/**/
			HitCollection *hitCollection;

			/**/
			TVectorD eigenValues;

			/**/
			TMatrixD eigenVectors;

			/*! Initialize the algorithm, i.e by initializing specific variables */
			virtual void Init();

			/*! Execute the algorithm */
			virtual void Execute();

			/*! Finalize the algorithm*/
			virtual void End();

			/*! Allow to check if everything is well set in the algorithm before starting it */
			virtual Return CheckConsistency();



	};  // class

}  // namespace 

#endif  //  PRINCIPALCOMPONENTANALYSIS_HH
