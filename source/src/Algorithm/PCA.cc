  /// \file PCA.cc
/*
 *
 * PCA.cc source template generated by fclass
 * Creation date : ven. juil. 12 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "Algorithm/PCA.hh"


using namespace std;

namespace baboon {

	PCA::PCA()
		: AbstractAlgorithm("PCA") {

		needData = false;
	}

	PCA::~PCA() {

		pcaRows.clear();
	}


	Return PCA::Init() {

		eigenValues.Clear();
		eigenVectors.Clear();

		return BABOON_SUCCESS();
	}


	Return PCA::CheckConsistency() {

		if( pcaRows.empty() )
			return BABOON_NOT_INITIALIZED("No rows in PCA");

		return BABOON_SUCCESS();
	}


	Return PCA::Execute() {


		DoubleVector rowMeans(pcaRows.size(),0);
		TMatrixD covMatrix(pcaRows.size(),pcaRows.size());

		for( unsigned int r=0 ; r<pcaRows.size() ; r++ )
			rowMeans.at(r) = accumulate ( pcaRows.at(r).begin() , pcaRows.at(r).end() , 0.0 ) *1.0 / pcaRows.at(r).size();



		for( unsigned int entry=0 ; entry<pcaRows.at(0).size() ; entry++ ) {

			for( unsigned int r1=0 ; r1<pcaRows.size() ; r1++ ) {

				for( unsigned int r2=r1 ; r2<pcaRows.size() ; r2++ ) {

					covMatrix( r1 , r2 ) += ( pcaRows.at(r1).at(entry) - rowMeans.at(r1) )*( pcaRows.at(r2).at(entry) - rowMeans.at(r2) ) / (pcaRows.at(0).size() - 1);

					// Covariance matrix is symmetric
					if( r1 != r2 )
						covMatrix( r2 , r1 ) = covMatrix( r1 , r2 );
				}
			}
		}

		TMatrixDEigen matrixEigen(covMatrix);

		eigenValues.ResizeTo( pcaRows.size() );
		eigenVectors.ResizeTo( pcaRows.size() , pcaRows.size() );
		eigenValues = matrixEigen.GetEigenValuesRe();
		eigenVectors = matrixEigen.GetEigenVectors();


		return BABOON_SUCCESS();
	}


	Return PCA::End() {

		pcaRows.clear();
		return BABOON_SUCCESS();
	}

	Return PCA::AddRow( const Row &row ) {

		if( row.empty() )
			return BABOON_INVALID_PARAMETER("PCA row is empty");

		if( pcaRows.empty() ) {

			pcaRows.push_back( row );
			return BABOON_SUCCESS();
		}
		else {
			if( row.size() != pcaRows.at(0).size() ) {
				return BABOON_INVALID_PARAMETER("PCA row has an incompatible size");
			}
			else {
				pcaRows.push_back( row );
				return BABOON_SUCCESS();
			}
		}
	}

}  // namespace 

