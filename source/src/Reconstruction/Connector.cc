  /// \file Connector.cc
/*
 *
 * Connector.cc source template generated by fclass
 * Creation date : jeu. oct. 10 2013
 *
 * This file is part of XXX libraries.
 * 
 * XXX is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * XXX is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with XXX.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author : Et� R�mi
 * @version
 * @copyright
 *
 *
 */


#include "Reconstruction/Connector.hh"

namespace baboon {

	template<typename T,typename S>
	Point<T,S>::Point() {

//		connectors = new ConnectorCollection< T , S >::type;
	}

	template<typename T,typename S>
	Point<T,S>::~Point() {

		// no connector deletion since, there
		// are not created in the Point object
		connectors.clear();
//		delete connectors;
	}

	template<typename T,typename S>
	void Point<T,S>::AddConnector( Connector<T,S> *connector ) {

		if( connector == 0 )
			return;

		if( !this->IsConnectedTo( connector ) ) {
			connectors.push_back( connector );
		}
		return;

	}

	template<typename T,typename S>
	void Point<T,S>::RemoveConnector( Connector<T,S> *connector ) {

		if( connector == 0 )
			return;

		typename ConnectorCollection< T , S >::type::iterator it = std::find( connectors.begin() , connectors.end() , connector );
		if( it != connectors.end() ) {
			connectors.erase( it );
		}
		return;

	}

	template<typename T,typename S>
	bool Point<T,S>::IsConnectedTo( Connector<T,S> *connector ) {

		return ( std::find( connectors.begin() , connectors.end() , connector ) != connectors.end() );
	}


	template<typename T,typename S>
	bool Point<T,S>::IsConnectedTo( Point<S,T> *point ) {

		for( unsigned int c=0 ; c<connectors.size() ; c++ )
			if( connectors.at( c )->GetSecond() == point )
				return true;

		return false;
	}


	template<typename T,typename S>
	void Point<T,S>::SetObject( T obj ) {

		object = obj;
	}

	template<typename T,typename S>
	T Point<T,S>::GetObject() {

		return object;
	}


	template<typename T,typename S>
	bool Point<T,S>::HasNoConnection() {

		return connectors.empty();
	}

	template<typename T,typename S>
	typename ConnectorCollection< T , S >::type &Point<T,S>::GetConnectors() {

		return connectors;
	}

// ----------------------------------------------------------------------------------------------------

	template<typename T , typename S>
	Connector<T,S>::Connector() {

		pointPair.first = 0;
		pointPair.second = 0;
	}

	template<typename T , typename S>
	Connector<T,S>::~Connector() {

		pointPair.first = 0;
		pointPair.second = 0;
	}

	template< typename T , typename S >
	void Connector<T,S>::Connect( Point<T,S> *point1 , Point<S,T> *point2 ) {

		if( point1 == 0 || point2 == 0 )
			return;

		pointPair.first = point1;
		pointPair.second = point2;
	}

	template< typename T , typename S >
	void Connector<T,S>::Disconnect() {

		pointPair.first = 0;
		pointPair.second = 0;
	}

	template<typename T,typename S>
	Point<T,S> *Connector<T,S>::GetFirst() {

		return pointPair.first;
	}

	template<typename T,typename S>
	Point<S,T> *Connector<T,S>::GetSecond() {

		return pointPair.second;
	}

	// for Point
	template baboon::Point< CaloHit* , CaloHit *>::Point();
	template baboon::Point< CaloHit *, CaloHit *>::~Point();
	template void Point< CaloHit * , CaloHit * >::AddConnector( Connector< CaloHit *, CaloHit* > *connector );
	template void baboon::Point< CaloHit * , CaloHit * >::RemoveConnector( baboon::Connector< CaloHit * , CaloHit * > *connector );
	template bool baboon::Point< CaloHit * , CaloHit * >::IsConnectedTo( baboon::Connector< CaloHit * , CaloHit * > *connector );
	template bool baboon::Point< CaloHit * , CaloHit * >::IsConnectedTo( Point< CaloHit * , CaloHit * > *point );
	template void baboon::Point< CaloHit * , CaloHit * >::SetObject( CaloHit *obj );
	template CaloHit *baboon::Point< CaloHit * , CaloHit * >::GetObject();
	template typename ConnectorCollection< CaloHit * , CaloHit * >::type &baboon::Point< CaloHit * , CaloHit * >::GetConnectors();

	// for Connector
	template void Connector< CaloHit * , CaloHit * >::Connect( Point< CaloHit * , CaloHit * > *point1 , Point< CaloHit * , CaloHit *> *point2 );
	template Connector< CaloHit * , CaloHit * >::Connector();
	template void Connector< CaloHit * , CaloHit * >::Disconnect();
	template Point< CaloHit * , CaloHit * > *Connector< CaloHit * , CaloHit * >::GetFirst();
	template Point< CaloHit * , CaloHit * > *Connector< CaloHit * , CaloHit * >::GetSecond();


}  // namespace 

