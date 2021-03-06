  /// \file Connector.hh
/*
 *
 * Connector.hh header template generated by fclass
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


#ifndef CONNECTOR_HH
#define CONNECTOR_HH

#include <iostream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>


namespace baboon {

	// forward declaration of Connector class
	template<typename T, typename S>
	class Connector;

	// forward declaration of Point class
	template<typename T, typename S>
	class Point;

	// forward declaration of ConnectorCollection typedef
	template<typename T,typename S>
	struct ConnectorCollection {
		 typedef std::vector< Connector<T,S> * > type;
	};

	// forward declaration of PointPair typedef
	template<typename T,typename S>
	struct PointPair {
		typedef std::pair< Point<T,S> * , Point<S,T> * > type;
	};

	// forward declaration of PointCollection typedef
	template<typename T,typename S>
	struct  PointCollection {
		typedef std::vector< class Point<T,S> * > type;
	};

	// forward declaration of OrderedPointCollection typedef
	template<typename T,typename S>
	struct OrderedPointCollection {
		typedef std::map< int , class PointCollection<T,S>::type * > type;
	};

	/*!
	 *
	 * @brief Point class
	 *
	 */
	template<typename T,typename S>
	class Point {

		public:

			/*!
			*
			* @brief  Default constructor
			*
			*/
			Point() {

			}

			/*!
			*
			* @brief  Default destructor
			*
			*/
			virtual ~Point() {

				connectors.clear();
			}

			/*!
			 *
			 *
			 *
			 */
			void AddConnector( Connector<T,S> *connector ) {

				if( connector == 0 )
					return;

				if( !this->IsConnectedTo( connector ) ) {
					connectors.push_back( connector );
				}
				return;
			}

			/*!
			 *
			 *
			 *
			 */
			void RemoveConnector( Connector<T,S> *connector ) {

				if( connector == 0 )
					return;

				typename ConnectorCollection< T , S >::type::iterator it = std::find( connectors.begin() , connectors.end() , connector );
				if( it != connectors.end() ) {
					connectors.erase( it );
				}
				return;
			}

			/*!
			 *
			 *
			 *
			 */
			bool IsConnectedTo( Connector<T,S> *connector ) {

				return ( std::find( connectors.begin() , connectors.end() , connector ) != connectors.end() );
			}

			/*!
			 *
			 * @brief Return true if the point is connected to 'point'
			 *
			 */
			bool IsConnectedTo( Point<S,T> *point ) {

				for( unsigned int c=0 ; c<connectors.size() ; c++ )
					if( connectors.at( c )->GetSecond() == point )
						return true;

				return false;
			}

			/*!
			 *
			 * Set the object
			 *
			 */
			void SetObject( T obj ) {

				object = obj;
			}

			/*!
			 *
			 * @brief Return the object
			 *
			 */
			T GetObject() {

				return object;
			}

			/*!
			 *
			 * @brief Return true if there is no connection
			 *
			 */
			bool HasNoConnection() {

				return connectors.empty();
			}

			/*!
			 *
			 *  @brief Return the connectors
			 *
			 */
			typename ConnectorCollection< T , S >::type &GetConnectors() {

				return connectors;
			}

		protected:

			typename ConnectorCollection< T , S >::type connectors;          ///< The connectors
			T object;                                     	                   ///< The associated object

	};

	/*!
	 *
	 * @brief Connector class.
	 *
	 */
	template<typename T,typename S>
	class Connector {

		public:

			/*!
			 *
			 * @brief Default constructor
			 *
			 */
			Connector() {

				pointPair.first = 0;
				pointPair.second = 0;
				weight = 1.;
				isGood = true;
			}

			/*!
			 *
			 * @brief Default destructor
			 *
			 */
			~Connector() {

				pointPair.first = 0;
				pointPair.second = 0;
				weight = 1.;
			}

			/*
			 *
			 * @brief Connects two points
			 *
			 */
			void Connect( Point<T,S> *point1 , Point<S,T> *point2 ) {

				if( point1 == 0 || point2 == 0 )
					return;

				pointPair.first = point1;
				pointPair.second = point2;
				isGood = true;
			}

			/*!
			 *
			 *
			 *
			 */
			void Connect( Point<T,S> *point1 , Point<S,T> *point2 , const double &w ) {

				if( point1 == 0 || point2 == 0 )
					return;

				pointPair.first = point1;
				pointPair.second = point2;
				weight = w;
				isGood = true;
			}

			/*!
			 *
			 * @brief Disconnect the two points
			 *
			 */
			void Disconnect() {

				pointPair.first = 0;
				pointPair.second = 0;
				weight = 1.;
				isGood = false;
			}

			/*!
			 *
			 * @brief Return the first connected point
			 *
			 */
			Point<T,S> *GetFirst() {

				return pointPair.first;
			}

			/*!
			 *
			 * @brief Return the second connected point
			 *
			 */
			Point<S,T> *GetSecond() {

				return pointPair.second;
			}

			/*!
			 *
			 * @brief Return the weight of the connection
			 *
			 */
			const double &GetWeight() {

				return weight;
			}

			/*!
			 *
			 * @brief Set the connector a 'good' one
			 *
			 */
			void SetGood( bool b ) {

				isGood = b;
			}

			/*!
			 *
			 * @brief True if the connector is a 'good' connector
			 *
			 */
			bool IsGood() {

				return isGood;
			}

		protected:

			typename PointPair<T,S>::type pointPair;    ///< The point pair that are connected
			double weight;
			bool isGood;

	}; // class

}  // namespace 

#endif  //  CONNECTOR_HH
