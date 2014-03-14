  /// \file ExternInfoManager.hh
/*
 *
 * ExternInfoManager.hh header template generated by fclass
 * Creation date : mer. f�vr. 5 2014
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
 * @author : Remi Ete
 * @version
 * @copyright
 *
 *
 */


#ifndef EXTERNINFOMANAGER_HH
#define EXTERNINFOMANAGER_HH

#include <vector>

#include "Objects/TrackInfo.hh"

namespace baboon {

	/**
		* @brief  ExternInfoManager class
		*/
	class ExternInfoManager {

	public:

		/**
			* @brief Constructor
			*/
		static ExternInfoManager *GetInstance();

		/**
			* @brief Destructor
			*/
		static void Kill();

		/**
		 *
		 */
		TrackInfo *CreateTrackInfo();

		/**
		 *
		 */
		void RemoveTrackInfo( TrackInfo * );

		/**
		 *
		 */
		std::vector< TrackInfo * > &GetTrackInfos();

		/**
		 *
		 */
		const std::vector< TrackInfo * > &GetTrackInfos() const;

		/**
		 *
		 */
		void ClearAllContent();

		/**
		 *
		 */
		void ClearTrackInfoContent();

	protected:

		static ExternInfoManager* instance;
		ExternInfoManager();
		virtual ~ExternInfoManager();

		std::vector< TrackInfo * > trackInfos;

	};  // class

}  // namespace 

#endif  //  EXTERNINFOMANAGER_HH