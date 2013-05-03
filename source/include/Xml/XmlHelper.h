/**
 *  @file   PandoraPFANew/Framework/include/Helpers/XmlHelper.h
 * 
 *  @brief  Header file for the xml helper class.
 * 
 *  $Log: $
 */
#ifndef XMLHELPER_H
#define XMLHELPER_H

#include "Utilities/Globals.hh"
#include "Utilities/ReturnValues.hh"
#include "Xml/tinyxml.h"



namespace baboon {

	/**
	 *  XmlHelper class.
	 */
	class XmlHelper {

		public:

			template <typename T>
			static Return ReadValue(const TiXmlHandle &xmlHandle, const std::string &xmlElementName, T &t);


			template <typename T>
			static Return ReadVectorOfValues(const TiXmlHandle &xmlHandle, const std::string &xmlElementName, std::vector<T> &vector);


			template <typename T>
			static Return Read2DVectorOfValues(const TiXmlHandle &xmlHandle, const std::string &xmlElementName, const std::string &rowName,
				std::vector< std::vector<T> > &vector);


			static void TokenizeString(const std::string &inputString, StringVector &tokens, const std::string &delimiter = " ");

	};



} // namespace

#endif  //  XMLHELPER_H
