/*
 *
 * Exception.hh header template generated by fclass
 * Creation date : Thu Mar 14 23:07:11 2013
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 *
 * Written by : R. Eté
 */


#ifndef SDHCAL_EXCEPTION_HH
#define SDHCAL_EXCEPTION_HH

#include <iostream>
#include <string>
#include <cstdlib>
#include <exception>

namespace baboon {



	/*!
	 * @brief : Class Exception. The base class for all exceptions in baboon.
	 */
	class Exception : public std::exception {

	protected :

		std::string message;


		/*!
		 * @brief : Default Constructor
		 */
		Exception() throw() {;}

	public :

		/*!
		 * @brief : Constructor thrown
		 */
		Exception( const std::string& msg ) {
			message = "baboon::Exception : " + msg;
		}


		/*!
		 * @brief : Default Destructor
		 */
		virtual ~Exception() throw() {;}


		/*!
		 * @brief : return the exception message
		 */
		virtual const char* what() const throw() { return message.c_str(); }

	};



	/*!
	 * @brief : Class AlgorithmException.
	 */
	class AlgorithmException : public Exception {

	protected :

		/*!
		 * @brief : Default Constructor
		 */
		AlgorithmException() throw() {;}

	public :

		/*!
		 * @brief : Constructor thrown
		 */
		AlgorithmException( const std::string& msg ) {
			message = "sdhcal::AlgorithmException : " + msg;
		}


		/*!
		 * @brief : Default Destructor
		 */
		virtual ~AlgorithmException() throw() {;}


		/*!
		 * @brief : return the exception message
		 */
		virtual const char* what() const throw() { return message.c_str(); }

	};





	/*!
	 * @brief : Class AnalysisException.
	 */
	class AnalysisException : public Exception {

	protected :

		/*!
		 * @brief : Default Constructor
		 */
		AnalysisException() throw() {;}

	public :

		/*!
		 * @brief : Constructor thrown
		 */
		AnalysisException( const std::string& msg ) {
			message = "sdhcal::AnalysisException : " + msg;
		}


		/*!
		 * @brief : Default Destructor
		 */
		virtual ~AnalysisException() throw() {;}


		/*!
		 * @brief : return the exception message
		 */
		virtual const char* what() const throw() { return message.c_str(); }

	};

}


#endif  // SDHCAL_EXCEPTION_HH
