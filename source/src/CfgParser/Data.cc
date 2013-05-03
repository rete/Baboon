/*
 *
 * //Data storage for configuration file parser
 *
 * Use and copying of this software and preparation of derivative works
 * based upon this software are permitted. Any copy of this software or
 * of any derivative work must include this paragraph.
 *
 * Written by : rete
 *
 */

#include "CfgParser/Data.hh"



using namespace std;



namespace cfgparser {


	Data::Data() {}


	Data::Data( const string &t)
		: type(t) {}


	Data::Data( Data const &data ) {

		type = data.type;
		keyValueMap = data.keyValueMap;
	}

	Data::~Data() {}


	void Data::Append( const string &k , const string &v ) {

		keyValueMap[k] = v;
	}


	bool Data::HasKey( const string &k ) {

		KeyValueMap::iterator it;
		for (it=keyValueMap.begin() ; it!=keyValueMap.end() ; it++)
			if(it->first==k) return true;
		return false;
	}


	void Data::GetValue( const string &k , string *v ) {

		if( k.empty() )
			throw DataException("Key empty while getting value!");
		string key = k;
		NormalizeName(&key);
		bool found = false;
		KeyValueMap::iterator it;

		if( keyValueMap.empty() ) return;
		for (it=keyValueMap.begin() ; it!=keyValueMap.end() ; it++) {
			if( it->first == key ) {
				*v = it->second;
				found = true;
			}
		}

		if(!found)
			throw DataException("Key \'"+key+"\' not found in data type \'"+this->GetType()+"\'");
	}


	void Data::GetValue( const string &k, int *v){

		string value;
		this->GetValue(k,&value);
		*v = atoi(value.c_str());
	}


	void Data::GetValue( const string &k , double *v){

		string value;
		this->GetValue(k,&value);
		*v = atof(value.c_str());
	}


	void Data::GetValue( const string &k , bool *v) {

		string value;
		this->GetValue(k,&value);

		if(value == "true" || value == "yes" || value == "1" || value == "on")
			*v = true;
		else if(value == "false" || value == "no" || value == "0" || value == "off")
			*v = false;
		else throw DataException("\'" + value + "\' is not a valid boolean value (key \'" + k + "\')");

	}


	void Data::GetValue( const string &k , vector<string> *v ){

		string value;
		this->GetValue(k,&value);
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



	void Data::GetValue( const std::string &k , std::vector<double> *v ){

		string value;
		GetValue( k , &value);
		string s;
		char motif = ':';
		for (unsigned int i=0 ; i<value.size() ; i++)
		{
			if(value[i] != motif) s.push_back(value[i]);
			else {v->push_back( atof( s.c_str() ) ); s = "";}
			if(i==value.size()-1) v->push_back( atof( s.c_str() ) );
		}
		if(v->empty()) throw DataException("Vector value empty (key \'"+k+"\')");
	}





	void Data::GetValue( const std::string &k , std::vector<int> *v ){

		string value;
		GetValue( k , &value);
		string s;
		char motif = ':';
		for (unsigned int i=0 ; i<value.size() ; i++)
		{
			if(value[i] != motif) s.push_back(value[i]);
			else {v->push_back( atoi( s.c_str() ) ); s = "";}
			if(i==value.size()-1) v->push_back( atoi( s.c_str() ) );
		}
		if(v->empty()) throw DataException("Vector value empty (key \'"+k+"\')");
	}




	bool Data::IsEmpty() {
		return keyValueMap.empty();
	}


	void Data::PrintData() {

		if(this->IsEmpty()) return;

		else {
			cout << "data is not empty !" << endl;
			cout << "[" << type << "]" << endl;
			KeyValueMap::iterator it;
			if(keyValueMap.size() != 0)
				for(it=keyValueMap.begin() ; it!=keyValueMap.end() ; it++)
					cout << "  " << (*it).first << " = " << (*it).second << endl;
		}
	}


	Data& Data::operator +=( const Data& data ) {

		KeyValueMap::iterator it;
		for (it = data.GetKeyValueMap().begin() ; it != data.GetKeyValueMap().end() ; it++)
			this->Append((*it).first,(*it).second);
		return *this;
	}


	Data operator+ ( const Data &data1 , const Data &data2 ) {

		Data data;
		KeyValueMap map1 = data1.GetKeyValueMap();
		KeyValueMap map2 = data2.GetKeyValueMap();
		KeyValueMap::iterator it;
		for (it = map1.begin() ; it != map1.end() ; it++)
			data.Append((*it).first,(*it).second);

		for (it = map2.begin() ; it != map2.end() ; it++)
			data.Append((*it).first,(*it).second);

		return data;
	}


	void Data::Clear() {
		keyValueMap.clear();
	}



}  //  end namespace cfgparser
