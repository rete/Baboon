  /// \file RawCfgParser.cc
/*
 *
 * RawCfgParser.cc source template generated by fclass
 * Creation date : jeu. juin 6 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#include "CfgParser/RawCfgParser.hh"

using namespace std;

namespace cfgparser {



	RawCfgParser::RawCfgParser() {

		sections.clear();
	}


	RawCfgParser::~RawCfgParser() {


		sections.clear();
	}


	StatusCode RawCfgParser::CreateSection( const std::string &secName ) {

		if( !this->HasSection( secName ) ) {
			Section *section = new Section( secName );
			sections.push_back( section );
//			cout << "success returned while creating a section" << endl;
//			cout << "nb of registered sections : " << sections.size() << endl;
			return CFGPARSER_SUCCESS();
		}
		return CFGPARSER_DUPLICATE_SECTION_ERROR( "Section " + secName + " already exists!" );
	}


	StatusCode RawCfgParser::AddSection( Section *section ) {

		if( section == 0 )
			return CFGPARSER_MISSING_SECTION_NAME_ERROR("Assertion section != 0 failed");

		if( !this->HasSection( section ) ) {
			sections.push_back( section );
			return CFGPARSER_SUCCESS();
		}
		else return CFGPARSER_DUPLICATE_SECTION_ERROR( "Section " + section->GetName() + " already exists!" );
	}


	StatusCode RawCfgParser::GetValue( const std::string& sectionName , const std::string &key, std::string *value ) const {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->GetValue( key , value ) );
		return CFGPARSER_SUCCESS();
	}


	StatusCode RawCfgParser::GetValue( const std::string& sectionName , const std::string &key, int *value ) const {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->GetValue( key , value ) );
		return CFGPARSER_SUCCESS();
	}


	StatusCode RawCfgParser::GetValue( const std::string& sectionName , const std::string &key, double *value ) const {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->GetValue( key , value ) );
		return CFGPARSER_SUCCESS();
	}


	StatusCode RawCfgParser::GetValue( const std::string& sectionName , const std::string &key, bool *value ) const {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->GetValue( key , value ) );
		return CFGPARSER_SUCCESS();
	}


	StatusCode RawCfgParser::GetValue( const std::string& sectionName , const std::string &key, std::vector< std::string > *value ) const {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->GetValue( key , value ) );
		return CFGPARSER_SUCCESS();
	}


	StatusCode RawCfgParser::GetValue( const std::string& sectionName , const std::string &key, std::vector< int > *value ) const {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->GetValue( key , value ) );
		return CFGPARSER_SUCCESS();
	}


	StatusCode RawCfgParser::GetValue( const std::string& sectionName , const std::string &key, std::vector< double > *value ) const {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->GetValue( key , value ) );
		return CFGPARSER_SUCCESS();
	}


	StatusCode RawCfgParser::GetValue( const std::string& sectionName , const std::string &key, std::vector< bool > *value ) const {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->GetValue( key , value ) );
		return CFGPARSER_SUCCESS();
	}


	bool RawCfgParser::RawCfgParser::HasKey( const std::string &sectionName , const std::string &key ) const {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		return section->HasKey( key );
	}


	bool RawCfgParser::HasSection( const std::string &secName ) const {

		for( unsigned int i=0 ; i<sections.size() ; i++ ) {
			if( sections.at(i)->GetName() == secName )
				return true;
		}
		return false;
	}


	bool RawCfgParser::HasSection( const Section *section ) const {

		for ( unsigned int i=0 ; i<sections.size() ; i++ ) {
			if( sections.at(i)->GetName() == section->GetName() )
				return true;
		}
		return false;
	}


	StatusCode RawCfgParser::GetSection( const std::string &sectionName , Section &sec ) const {

		if( sectionName.empty() )
			return CFGPARSER_MISSING_SECTION_NAME_ERROR("Assertion !sectionName.empty() failed");

		for( unsigned int i=0 ; i<sections.size() ; i++ ) {
			if( sections.at(i)->GetName() == sectionName ) {
				sec = *sections.at(i);
				return CFGPARSER_SUCCESS();
			}
		}
		return CFGPARSER_NO_SECTION_ERROR( "Section " + sectionName + " not found!" );
	}

	StatusCode RawCfgParser::GetSection( const std::string &sectionName , Section *&sec ) const {

		if( sectionName.empty() )
			return CFGPARSER_MISSING_SECTION_NAME_ERROR("Assertion !section.empty() failed");

		for( unsigned int i=0 ; i<sections.size() ; i++ ) {
			if( sections.at(i)->GetName() == sectionName ) {
				sec = sections.at(i);
				return CFGPARSER_SUCCESS();
			}
		}
		return CFGPARSER_NO_SECTION_ERROR( "Section " + sectionName + " not found!" );
	}

	StatusCode RawCfgParser::GetKeys( const std::string &sectionName , std::vector< std::string > &keys ) const {

		if( sectionName.empty() )
			return CFGPARSER_MISSING_SECTION_NAME_ERROR("Assertion !section.empty() failed");

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		KeyValueMap keyValueMap = section->GetKeyValueMap();
		keys.clear();
		KeyValueMap::iterator it = keyValueMap.begin();
		for( ; it!=keyValueMap.end() ; it++ )
			keys.push_back( it->first );

		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::Read( const std::string &fileName ) {

		if( fileName.empty() )
			return CFGPARSER_ERROR("Assertion !fileName.empty() failed");

		ifstream *cfgFile = new ifstream();
		cfgFile->open( fileName.c_str() );

		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , Read( cfgFile ) );

		cfgFile->close();
		delete cfgFile;

		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::Read( const std::vector< std::string > &fileNames ) {

		if( fileNames.empty() )
			return CFGPARSER_ERROR("Assertion !fileNames.empty() failed");

		for( unsigned int f=0 ; f<fileNames.size() ; f++ ) {
			CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->Read( fileNames.at( f ) ) );
		}

		return CFGPARSER_SUCCESS();
	}




	StatusCode RawCfgParser::Read( std::ifstream& stream ) {

		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->Read( &stream ) );
		return CFGPARSER_SUCCESS();
	}




	StatusCode RawCfgParser::Read( std::ifstream* stream ) {

		if( stream == 0 )
			return CFGPARSER_ERROR("Assertion stream != 0 failed");

		if( !stream->is_open() )
			return CFGPARSER_ERROR("Cfg file stream is not opened or doesn't exists!");

		char quoteOpen = 0;
		Section *currentSection = 0;


		while( !stream->eof() ) {

			string sectionName;
			string key;
			string value;
			string lineBuffer;
			char ch;
			stream->get( ch );

			lineBuffer.erase();
			while ( !stream->eof() && ch != '\n') {
				lineBuffer += ch;
				stream->get( ch );
			}

			if( lineBuffer.empty() )
				continue;

			// skip blank char
			int i = 0;
			for( i=0 ; i<lineBuffer.size() ; i++ )
				if( lineBuffer.at(i) != ' ' )
					break;

			// skip blank line
			if( lineBuffer.substr( i ).empty() )
				continue;

			lineBuffer = lineBuffer.substr( i );

			// skip comment lines
			if( lineBuffer.at(0) == '#' )
				continue;

			// beginning of a section
			if( lineBuffer.at(0) == '[' ) {

				// read the section name
				int j=1;
				sectionName.clear();
				while( lineBuffer.at( j ) != ']' ) {

					if( lineBuffer.at( j ) == ' ' ) {
						j++;
						continue;
					}
					sectionName += lineBuffer.at( j );
					j++;
				}

				// create it and get it as current section
				currentSection = new Section( sectionName );
				CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->AddSection( currentSection ) );
				continue;
			}

			if( currentSection == 0 )
				continue;

			// reading a key
			key.clear();
			for( i=0 ; i<lineBuffer.at(i) ; i++ ) {

				if( lineBuffer.at(i) == '=' )
					break;
				key += lineBuffer.at(i);
			}
			NormalizeName( &key );

			lineBuffer.erase(0,1);

			//reading a value
			value.clear();
			while (i < lineBuffer.length()) {

				if (quoteOpen == 0) {
					if (lineBuffer[i] == '#' || lineBuffer[i] == ';') break;
					if (lineBuffer[i] == '\'' || lineBuffer[i] == '"') quoteOpen = lineBuffer[i];
					if (lineBuffer[i] == '\\'
					&& i < lineBuffer.length()+1
					&& (lineBuffer[i+1] == '#' || lineBuffer[i+1] == ';') ) i++;
				}
				else if (lineBuffer[i] == quoteOpen) quoteOpen = 0;

				value += lineBuffer[i++];
			}
			StrTrim( &value );
			currentSection->SetValue( key , value );
		}
		return CFGPARSER_SUCCESS();
	}







	StatusCode RawCfgParser::RemoveKey( const std::string &sectionName , const std::string &key ) {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->RemoveKey( key ) );
		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::RemoveSection( const std::string &sectionName ) {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		sections.erase( std::find( sections.begin() , sections.end() , section ) );
		delete section;
		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::GetSections( SectionCollection &sec ) const {

		sec = this->sections;
		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::SetValue( const std::string& sectionName , const std::string &key, const std::string &value ) {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->SetValue( key , value ) );
		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::SetValue( const std::string& sectionName , const std::string &key, const int &value ) {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		ostringstream ss;
		ss << value;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->SetValue( key , ss.str() ) );
		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::SetValue( const std::string& sectionName , const std::string &key, const double &value ) {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		ostringstream ss;
		ss << value;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->SetValue( key , ss.str() ) );
		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::SetValue( const std::string& sectionName , const std::string &key, const bool &value ) {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );

		if( value ) {
			CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->SetValue( key , "true" ) );
		}
		else {
			CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->SetValue( key , "false" ) );
		}

		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::SetValue( const std::string& sectionName , const std::string &key, const std::vector< std::string > &value ) {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		ostringstream ss;
		for( unsigned int i=0 ; i<value.size() ; i++ ) {
			ss << value.at(i);
			if( i == value.size() - 1 )
				break;
			ss << ":";
		}
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->SetValue( key , ss.str() ) );
		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::SetValue( const std::string& sectionName , const std::string &key, const std::vector< int > &value ) {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		ostringstream ss;
		for( unsigned int i=0 ; i<value.size() ; i++ ) {
			ss << value.at(i);
			if( i == value.size() - 1 )
				break;
			ss << ":";
		}
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->SetValue( key , ss.str() ) );
		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::SetValue( const std::string& sectionName , const std::string &key, const std::vector< double > &value ) {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		ostringstream ss;
		for( unsigned int i=0 ; i<value.size() ; i++ ) {
			ss << value.at(i);
			if( i == value.size() - 1 )
				break;
			ss << ":";
		}
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->SetValue( key , ss.str() ) );
		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::SetValue( const std::string& sectionName , const std::string &key, const std::vector< bool > &value ) {

		Section *section = 0;
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->GetSection( sectionName , section ) );
		ostringstream ss;
		for( unsigned int i=0 ; i<value.size() ; i++ ) {
			if( value.at(i) )
				ss << "true";
			else ss << "false";
			if( i == value.size() - 1 )
				break;
			ss << ":";
		}
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , section->SetValue( key , ss.str() ) );
		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::Write( const std::string &fileName ) const {

		ofstream *cfgFile = new ofstream();
		cfgFile->open( fileName.c_str() );
		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->Write( cfgFile ) );
		cfgFile->close();
		delete cfgFile;
		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::Write( std::ofstream& stream ) const {

		CFGPARSER_THROW_RESULT_IF( CFGPARSER_SUCCESS() , != , this->Write( &stream ) );
		return CFGPARSER_SUCCESS();
	}

	StatusCode RawCfgParser::Write( std::ofstream* stream ) const {

		if( stream == 0 )
			return CFGPARSER_ERROR("Assertion stream != 0 failed");

		if( !stream->is_open() )
			return CFGPARSER_ERROR("Cfg file stream is not opened or doesn't exists!");

		if( sections.empty() )
			return CFGPARSER_ERROR("No section to write out!");

		time_t now = time( 0 );
		*stream << "######################################################" << endl;
		*stream << "# CfgParser libraries for configuration file parsing" << endl;
		*stream << "# Generated by a RawCfgParser" << endl;
		*stream << "# @date : " << ctime( &now );
		*stream << "######################################################" << endl;
		*stream << endl;
		*stream << endl;
		*stream << endl;

		for( unsigned int i=0 ; i<sections.size() ; i++ ) {

			*stream << "[" << sections.at(i)->GetName() << "]" << endl;
			*stream << endl;
			KeyValueMap keyValueMap = sections.at(i)->GetKeyValueMap();
			KeyValueMap::iterator it;
			for( it=keyValueMap.begin() ; it!=keyValueMap.end() ; it++ )
				*stream << it->first << " = " << it->second << endl;
			*stream << endl;
			*stream << endl;
		}
		*stream << endl;

		return CFGPARSER_SUCCESS();
	}


	StatusCode RawCfgParser::Print() const {

		if( sections.empty() )
			return CFGPARSER_SUCCESS();

		for( unsigned int i=0 ; i<sections.size() ; i++ ) {

			cout << "[" << sections.at(i)->GetName() << "]" << endl;
			cout << endl;
			KeyValueMap keyValueMap = sections.at(i)->GetKeyValueMap();
			KeyValueMap::iterator it;
			for( it=keyValueMap.begin() ; it!=keyValueMap.end() ; it++ )
				cout << it->first << " = " << it->second << endl;
			cout << endl;
		}
		return CFGPARSER_SUCCESS();
	}

}  // namespace 

