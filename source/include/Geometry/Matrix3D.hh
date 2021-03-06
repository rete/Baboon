/*
 *
 * Matrix3D.hh header template generated by fclass
 * Creation date : Thu Mar  7 16:46:22 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#ifndef MATRIX3D_HH
#define MATRIX3D_HH

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>


namespace baboon {

	/*!
	 * Class Matrix3D.
	 */
	template <typename T>
	class Matrix3D {

		protected :

			/*! buffer containing all values */
			T *buffer;
			/*! x size of the 3D matrix */
			unsigned int sizeX;
			/*! y size of the 3D matrix */
			unsigned int sizeY;
			/*! z size of the 3D matrix */
			unsigned int sizeZ;

		public :

			/*! Default Constructor. Initialize all values to 0 */
			Matrix3D() { buffer = NULL; }

			void Initialize( T* p , unsigned int x_size , unsigned int y_size , unsigned int z_size ) {
				buffer = p;
				sizeX = x_size;
				sizeY = y_size;
				sizeZ = z_size;
				this->Clear();
			}

			/*! Default Destructor */
			virtual ~Matrix3D()
				{ if(buffer!=0) delete buffer; }

			/*! */
			void Fill(T* p)
				{ memcpy( buffer, p, sizeX*sizeY*sizeZ*sizeof( T ) ); }

			/*! clear all the element in the buffer */
			void Clear()
				{ memset( buffer, 0, sizeX*sizeY*sizeZ*sizeof( T ) ); }

			/*! return the buffer */
			inline T* GetPtr() { return buffer; }

			/*! return the value at (i,j,k) */
			T  GetValue(unsigned int i, unsigned int j, unsigned int k)
				{ return buffer[ ( i*sizeY + j )*sizeZ + k ]; }

			/*! set the value at (i,j,k) */
			void SetValue(unsigned int i, unsigned int j, unsigned int k, T val)
				{ buffer[ (i*sizeY + j)*sizeZ + k ] = val; }

			/*! return the x size of the matrix */
			inline unsigned int GetXSize()
				{ return sizeX; }

			/*! return the y size of the matrix */
			inline unsigned int GetYSize()
				{ return sizeY; }

			/*! return the z size of the matrix */
			inline unsigned int GetZSize()
				{ return sizeZ; }

			int HowManyOf( T val ) {
				int count = 0;
				for( unsigned int i=0 ; i<this->GetXSize() ; i++ ) {
					for( unsigned int j=0 ; j<this->GetYSize() ; j++ ) {
						for( unsigned int k=0 ; k<this->GetZSize() ; k++ ) {
							if( this->GetValue(i,j,k) == val ) count++;
						}
					}
				}
				return count;
			}

	};

}

#endif  // MATRIX3D_HH
