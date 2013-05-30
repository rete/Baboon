/**
 *  @file   PandoraPFANew/Framework/src/Helpers/XmlHelper.cc
 * 
 *  @brief  Implementation of the xml helper class.
 * 
 *  $Log: $
 */


#include "Xml/XmlHelper.h"



namespace baboon {

	template <typename T>
	inline Return XmlHelper::ReadValue(const TiXmlHandle &xmlHandle, const std::string &xmlElementName, T &t) {

		const TiXmlElement *const pXmlElement = xmlHandle.FirstChild(xmlElementName).Element();

		if (NULL == pXmlElement)
			return BABOON_NOT_INITIALIZED("NULL == pXmlElement");

		if (!StringToType(pXmlElement->GetText(), t))
			return BABOON_ERROR("Can't convert string to asked type");

		return BABOON_SUCCESS();
	}

	template <>
	inline Return XmlHelper::ReadValue<bool>(const TiXmlHandle &xmlHandle, const std::string &xmlElementName, bool &t) {

		const TiXmlElement *const pXmlElement = xmlHandle.FirstChild(xmlElementName).Element();

		if (NULL == pXmlElement)
			return BABOON_NOT_INITIALIZED("NULL == pXmlElement");

		const std::string xmlElementString = pXmlElement->GetText();

		if ((xmlElementString == "1") || (xmlElementString == "true")) {

			t = true;
		}
		else if ((xmlElementString == "0") || (xmlElementString == "false")) {

			t = false;
		}
		else {
			return BABOON_ERROR("Invalid boolean value");
		}

		return BABOON_SUCCESS();
	}

	//------------------------------------------------------------------------------------------------------------------------------------------

	template <typename T>
	inline Return XmlHelper::ReadVectorOfValues(const TiXmlHandle &xmlHandle, const std::string &xmlElementName, std::vector<T> &vector) {

		const TiXmlElement *const pXmlElement = xmlHandle.FirstChild(xmlElementName).Element();

		if (NULL == pXmlElement)
			return BABOON_NOT_INITIALIZED("NULL == pXmlElement");

		StringVector tokens;
		TokenizeString(pXmlElement->GetText(), tokens);

		for (StringVector::const_iterator iter = tokens.begin(), iterEnd = tokens.end(); iter != iterEnd; ++iter)
		{
			T t;

			if (!StringToType(*iter, t))
				return BABOON_ERROR("Can't convert string to asked type");

			vector.push_back(t);
		}

		return BABOON_SUCCESS();
	}

	//------------------------------------------------------------------------------------------------------------------------------------------

	template <typename T>
	inline Return XmlHelper::Read2DVectorOfValues(const TiXmlHandle &xmlHandle, const std::string &xmlElementName, const std::string &rowName,
		std::vector< std::vector<T> > &vector)
	{
		TiXmlElement *pXmlElement = xmlHandle.FirstChild(xmlElementName).Element();

		if (NULL == pXmlElement)
			return BABOON_NOT_INITIALIZED("NULL == pXmlElement");

		TiXmlElement *pXmlRowElement = TiXmlHandle(pXmlElement).FirstChild(rowName).Element();

		if (NULL == pXmlRowElement)
			return BABOON_NOT_INITIALIZED("NULL == pXmlRowElement");

		for ( ; NULL != pXmlRowElement; pXmlRowElement = pXmlRowElement->NextSiblingElement(rowName))
		{
			std::vector<T> rowVector;

			StringVector tokens;
			TokenizeString(pXmlRowElement->GetText(), tokens);

			for (StringVector::const_iterator iter = tokens.begin(), iterEnd = tokens.end(); iter != iterEnd; ++iter)
			{
				T t;

				if (!StringToType(*iter, t))
					return BABOON_ERROR("Can't convert string to asked type");

				rowVector.push_back(t);
			}

			vector.push_back(rowVector);
		}

		return BABOON_SUCCESS();
	}

	void XmlHelper::TokenizeString(const std::string &inputString, StringVector &tokens, const std::string &delimiter) {

		std::string::size_type lastPos = inputString.find_first_not_of(delimiter, 0);
		std::string::size_type pos     = inputString.find_first_of(delimiter, lastPos);

		while ((std::string::npos != pos) || (std::string::npos != lastPos)) {

			tokens.push_back(inputString.substr(lastPos, pos - lastPos));
			lastPos = inputString.find_first_not_of(delimiter, pos);
			pos = inputString.find_first_of(delimiter, lastPos);
		}
	}

} // namespace pandora
