/// \file CfgParser.cc
/*
 *
 * @brief : Generic configuration file parser
 *
 * Use and copying of this software and preparation of derivative works
 * based upon this software are permitted. Any copy of this software or
 * of any derivative work must include this paragraph.
 *
 * @author : rete
 *
 */


using namespace std;


#include "CfgParser/CfgParser.hh"


//////////////////////////////////////////////////////////////////////////
////////////// DON'T TOUCH IF DON'T KNOW WHAT YOU ARE DOING //////////////
//////////////////////////////////////////////////////////////////////////
//
//                     / \  //\
//      |\___/|      /   \//  \\
//      /0  0  \__  /    //  | \ \
//     /     /  \/_/    //   |  \  \
//     @_^_@'/   \/_   //    |   \   \
//     //_^_/     \/_ //     |    \    \
//  ( //) |        \///      |     \     \
//( / /) _|_ /   )  //       |      \     _\
//( // /) '/,_ _ _/  ( ; -.    |    _ _\.-~        .-~~~^-.
//(( / / )) ,-{        _      `-.|.-~-.           .~         `.
//(( // / ))  '/\      /                 ~-. _ .-~      .-~^-.  \
//(( /// ))      `.   {            }                   /      \  \
//(( / ))     .----~-.\        \-'                 .~         \  `. \^-.
//       ///.----..>        \             _ -~             `.  ^-`  ^-_
//         ///-._ _ _ _ _ _ _}^ - - - - ~                     ~-- ,.-~
//                                                            /.-~


namespace cfgparser {



	CfgParser::CfgParser() {

		configFileName = "";
		file = new fstream();
		alreadyRead = false;
	}

	CfgParser::CfgParser( const string &fileName) {

		configFileName = fileName;
		file = new fstream();
		alreadyRead = false;
	}


	CfgParser::~CfgParser() {

		dataList.clear();
		delete file;
	}



	void CfgParser::Read()
	throw ( CfgParserException , std::exception ) {

		if( alreadyRead ) return;

		if( configFileName.empty() )
			throw CfgParserException( "While reading : cfg file name is empty."
					" Please set it before read it!" );

		file->open( configFileName.c_str() , ios::in );

		if( file->fail() )
			throw CfgParserException("\'"+configFileName+"\' configuration file doesn't exists or can't be opened!");

		dataList.clear();

		char quoteOpen = 0;
		string buf, lastType;
		Data data;

		// Read the whole file
		while(!file->eof() ) {
			int ch = file->get();

			string type;
			string key, value;

			// Get a line in the file.
			buf.erase();
			while (!file->eof() && ch != '\n') {
				buf += ch;
				ch = file->get();
			}

			if(file->eof()) {
				dataList.push_back(data);
				break;
			}
			unsigned int i = 0;
			while (i < buf.length() && isspace(buf[i])) ++i;
				/* read section */
			if (i < buf.length() && buf[i] == '[') {
				if(!data.IsEmpty()) dataList.push_back(data);
				data.Clear();
				/* section name */
				++i;
				while (i < buf.length() && buf[i] != ']')
					type += buf[i++];
				/* notice actual section */
				NormalizeName(&type);
				lastType = type;

				/* skip ']' */
				if (i < buf.length() && buf[i] == ']') ++i;
				/* skip intermediate spaces */
				while (i < buf.length() && isspace(buf[i])) ++i;
			}
			else type = lastType;

			bool containsSpace = false;
			for(unsigned int j=0 ; j<type.size() ; j++) if(type[j]==' ') containsSpace = true;

			if(containsSpace) throw CfgParserException("Data type '"+type+"' contains space(s). Data type mustn't contain space. Try something else.");

			data.SetType(type);

				/* entry name */
			while (i < buf.length() && buf[i] != '=' && buf[i] != '#'&& buf[i] != ';') key += buf[i++];

			NormalizeName(&key);
				/* regular content */
				/* equal sign */
			if (i < buf.length() && buf[i] == '=') {
				++i;
				/* read value */
				while (i < buf.length()) {
					/* are we inside quotes? */
					if (quoteOpen == 0) {
						if (buf[i] == '#' || buf[i] == ';') break;
						if (buf[i] == '\'' || buf[i] == '"') quoteOpen = buf[i];
						if (buf[i] == '\\'
						&& i < buf.length()+1
						&& (buf[i+1] == '#' || buf[i+1] == ';') ) i++;
					}
					else if (buf[i] == quoteOpen) quoteOpen = 0;

					value += buf[i++];
				}
				StrTrim(&value);
				data.Append(key,value);
				}
		} // End of reading

		file->close();
		alreadyRead = true;
	}





	Data CfgParser::GetData( const string &t ) {

		DataList::iterator it;
		bool typeFound = false;
		string type = t;
		NormalizeName(&type);

		if(dataList.size() != 0) {
			for(it=dataList.begin() ; it!=dataList.end() ; it++) {
				if( (it->GetType()==type)) {
					typeFound = true;
					if((*it).GetKeyValueMap().size() != 0) {
						return (*it);
					}
				}
			}
			if (!typeFound) throw CfgParserException("Data type '" + type + "' doesn't exists!");
		}
		else if (!alreadyRead) throw CfgParserException("File not read and data list is empty!");
		else throw CfgParserException("Data list empty! Try to read the file again");

		Data emptyDataObject;
		return emptyDataObject;
	}




	void CfgParser::PrintDataList() {

		if(dataList.size() != 0) {
			for(unsigned int i=0 ; i<dataList.size() ; i++) {
				dataList.at(i).PrintData();
				cout << endl;
			}
		}
		else if (!alreadyRead) throw CfgParserException("File not read and data list is empty!");
		else throw CfgParserException("Data list empty! Try to read the file again");
	}


	void CfgParser::GetValue( const string &t , const string &k , string *v ) {

		string key = k;
		NormalizeName(&key);
		DataList::iterator it;
		bool typeFound = false;

		if(dataList.size() != 0) {
			for(it=dataList.begin() ; it!=dataList.end() ; it++) {
				if(it->GetType() == t) {
					typeFound = true;
					it->GetValue(key,v);
				}
			}
			if (!typeFound) throw CfgParserException("Data type '" + t + "' doesn't exists!");
		}
		else if (!alreadyRead) throw CfgParserException("File not read and data list is empty!");
		else throw CfgParserException("Data list empty! Try to read the file again");
	}




	void CfgParser::GetValue(const string &t , const string &k , bool *v ) {

		string value;
		GetValue( t , k , &value );

		if(value == "true" || value == "on" || value == "1" || value == "yes") {
			*v = true;
		}
		else if(value == "false" || value == "off" || value == "0" || value == "no") {
			*v = false;
		}
		else throw CfgParserException("\'" + value + "\' is not a valid boolean value (key \'" + k + "\')");
	}



	void CfgParser::GetValue( const string &t , const string &k , int *v ) {

		string value;
		GetValue( t , k , &value );
		*v = atoi(value.c_str());
	}




	void CfgParser::GetValue( const string &t , const string &k , double *v ) {

		string value;
		GetValue( t , k , &value );
		*v = atof(value.c_str());
	}




	void CfgParser::GetValue( const std::string &t , const std::string &k , std::vector<std::string> *v ) {

		string value;
		GetValue( t , k , &value);
		string s;
		char motif = ':';
		for (unsigned int i=0 ; i<value.size() ; i++)
		{
			if(value[i] != motif) s.push_back(value[i]);
			else {v->push_back(s); s = "";}
			if(i==value.size()-1) v->push_back(s);
		}
		if(v->empty()) throw CfgParserException("Vector value empty (key \'"+k+"\')");
	}



	void CfgParser::GetValue( const std::string &t , const std::string &k , std::vector<double> *v ) {

		string value;
		GetValue( t , k , &value);
		string s;
		char motif = ':';
		for (unsigned int i=0 ; i<value.size() ; i++)
		{
			if(value[i] != motif) s.push_back(value[i]);
			else {v->push_back( atof( s.c_str() ) ); s = "";}
			if(i==value.size()-1) v->push_back( atof( s.c_str() ) );
		}
		if(v->empty()) throw CfgParserException("Vector value empty (key \'"+k+"\')");
	}



	void CfgParser::GetValue(const std::string &t , const std::string &k , std::vector<int> *v ) {

		string value;
		GetValue( t , k , &value);
		string s;
		char motif = ':';
		for (unsigned int i=0 ; i<value.size() ; i++)
		{
			if(value[i] != motif) s.push_back(value[i]);
			else {v->push_back( atoi( s.c_str() ) ); s = "";}
			if(i==value.size()-1) v->push_back( atoi( s.c_str() ) );
		}
		if(v->empty()) throw CfgParserException("Vector value empty (key \'"+k+"\')");
	}



}  //  end namespace cfgparser
