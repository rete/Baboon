  /// \file InputTree.hh
/*
 *
 * InputTree.hh header template generated by fclass
 * Creation date : dim. avr. 28 2013
 * Copyright (c) CNRS , IPNL
 *
 * All Right Reserved.
 * Use and copying of these libraries and preparation of derivative works
 * based upon these libraries are permitted. Any copy of these libraries
 * must include this copyright notice.
 * 
 * @author : rete
 */


#ifndef INPUTTREEWRAPPER_HH
#define INPUTTREEWRAPPER_HH


#include <iostream> 
#include <string> 
#include <cstdlib> 
#include <cmath> 
#include <vector> 
#include <utility>
#include <map>
#include <stdexcept>

// root includes
#include <TObjArray.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>


namespace baboon {


	typedef std::map< std::string , TBranch* > TBranchMap;

	/*!
	 *
	 * @brief : Class InputTreeWrapper
	 *
	 */

	class InputTTreeWrapper {


		private :

			// the input tree
			TTree *tree;

		public:

			/*!
			 *
			 * @ brief : Constructor
			 *
			 */
			InputTTreeWrapper( TTree *pTree );


			/*!
			 *
			 * @ brief : Destructor
			 *
			 */
			~InputTTreeWrapper();


			/*!
			 *
			 * @ brief : a simple print of the tree
			 *
			 */
			inline void Print() { tree->Print(); }


			/*!
			 *
			 * @ brief : a simple scan of the tree.
			 *
			 */
			inline void Scan() { tree->Scan(); }


			/*!
			 *
			 * @ brief : return tree tree instance
			 *
			 */
			inline TTree *GetTree() { return tree; }


			/*!
			 *
			 * @ brief : return the number of entries in the tree.
			 *
			 */
			inline int GetNbOfEntries() { return nbOfEntries; }


			/*!
			 *
			 * @ brief : Grab the current value of a given branch. Should ba called after LoadEntry(i);
			 *
			 */
			template <typename T>
			void GetValue( const std::string &branchName , T &val );

			/*!
			 *
			 * @ brief : Load an entry of the tree. Method to be called in a loop, at the begining.
			 *
			 */
			int LoadEntry( int entry );


		private :

			/*!
			 *
			 * @ brief : initialize the tree map
			 *
			 */
			void Init();


			/*!
			 *
			 * @ brief : Load the tree at entry
			 *
			 */
			int LoadTree( int entry );


			/*!
			 *
			 * @ brief : Get the entry nb 'entry'
			 *
			 */
			int GetEntry( int entry );

			TBranchMap *branchMap;
			int currentEntry;
			int nbOfEntries;

	};



}  // namespace 

#endif  //  INPUTTREEWRAPPER_HH
