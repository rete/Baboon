  /// \file Internal.hh
/*
 *
 * Internal.hh header template generated by fclass
 * Creation date : sam. mai 25 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef INTERNAL_HH
#define INTERNAL_HH

#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 

#include "Exception.hh"


#define BABOON_THROW_RESULT_IF(ReturnValue ,Operator ,FunctionResult )                                     \
	{                        																		        \
                                                     	 	 	 	 	 	 	 	 	 	 	 	 	 	\
		if( ReturnValue.fStatus Operator FunctionResult.fStatus ) { 										\
																											\
	        std::cerr << #FunctionResult << " throw " << FunctionResult.ToString() << std::endl;            \
	        std::cerr << "    in function : " << __FUNCTION__ << std::endl;                                 \
	        std::cerr << "    in file :     " << __FILE__ << " line#: " << __LINE__ << std::endl;           \
	        std::cerr << "    message :     " << FunctionResult.message << std::endl;                       \
	        throw ReturnValueException(FunctionResult.message);												\
		}																									\
	}


#define BABOON_THROW_RESULT( FunctionResult )                                                              \
	{																										\
		std::cerr << #FunctionResult << " throw " << FunctionResult.ToString() << std::endl;				\
		std::cerr << "    in function : " << __FUNCTION__ << std::endl;										\
		std::cerr << "    in file :     " << __FILE__ << std::endl;										\
		std::cerr << "    message :     " << FunctionResult.message << std::endl;										\
		throw ReturnValueException( FunctionResult.message );												\
	}


namespace baboon {

	/*!
	 * @brief : Class ReturnValueException.
	 */

	class ReturnValueException : public Exception {

		protected :

			/*!
			 * @brief : Default Constructor
			 */
			ReturnValueException() throw() {;}

		public :

			/*!
			 * @brief : Constructor thrown
			 */
			ReturnValueException( const std::string& msg ) {
				message = "baboon::ReturnValueException : " + msg;
			}


			/*!
			 * @brief : Default Destructor
			 */
			virtual ~ReturnValueException() throw() {;}


			/*!
			 * @brief : return the exception message
			 */
			virtual const char* what() const throw() { return message.c_str(); }

	};




}  // namespace 

#endif  //  INTERNAL_HH
