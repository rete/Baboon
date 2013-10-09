/*
 *
 * //LcioSplitFile.cc main file
 * Copyright (c) CNRS / IPNL
 * All Right Reserved.
 *
 * Use and copying of this software and preparation of derivative works
 * based upon this software are permitted. Any copy of this software or
 * of any derivative work must include this paragraph.
 *
 * Written by : R. Eté
 *
 *
 *
 *
 *
 *
 * Split a lcio file into N files.
 * The file is split according to the totalNbOfEvents/N.
 * If totalNbOfEvents/N is not an integer, the last file
 * will contains the remaining events. The files are renamed
 * with an additionnal suffixe "*_Ii.slcio".
 *
 *
 * Example :
 *    Split a file myFile.slcio with 10 events into 3 files
 *    leads to :
 *    - file 1 : 3 evts (myFile_I0.slcio)
 *    - file 2 : 3 evts (myFile_I1.slcio)
 *    - file 3 : 4 evts (myFile_I2.slcio)
 *
 *
 */


// std includes
#include <iostream>
#include <cstdlib>
#include <sstream>

// lcio includes
#include "IO/LCReader.h"
#include "IO/LCWriter.h"
#include "IOIMPL/LCFactory.h"
#include "EVENT/LCIO.h"
#include "lcio.h"

// tclap includes
#include "tclap/CmdLine.h"

using namespace std;
using namespace EVENT;
using namespace lcio;


int main ( int argc , char *argv[] ) {

	/*********************************
	 * Define the command line parser
	 *********************************/

	string footer = "Please report bug to <rete@ipnl.in2p3.fr>";
	TCLAP::CmdLine cmd(footer, ' ', "1.0");
	TCLAP::ValueArg<std::string> lcioFileArg(
				"f"
				,"slcio-file"
				,"lcio file to split"
				,true
				,""
				,"string" );
	TCLAP::ValueArg<unsigned int> numberOfFilesArg(
				"n"
				,"number-of-files"
				,"Number of of output files"
				,true
				,1
				,"unsigned int");

	cmd.add( lcioFileArg );
	cmd.add( numberOfFilesArg );
	cmd.parse( argc, argv );

	string lcioFile = lcioFileArg.getValue();
	unsigned int numberOfFiles = numberOfFilesArg.getValue();
	string fileRoot = lcioFile.substr( 0 , lcioFile.size() - 6 );
	string outputFileName = fileRoot + "_I0.slcio";

	IO::LCReader* lcReader = IOIMPL::LCFactory::getInstance()->createLCReader();
	IO::LCWriter* lcWriter = IOIMPL::LCFactory::getInstance()->createLCWriter();
	lcReader->open( lcioFile );
	lcWriter->open( outputFileName , LCIO::WRITE_NEW );
	EVENT::LCEvent *evt;

	int nbOfEvents = 0;
	int nbOfWroteFiles = 1;
	int nbOfWroteEvents = 0;

	while( (evt = lcReader->readNextEvent()) != 0  ) nbOfEvents++;

	cout << "File "<< lcioFile << " will be split in " << numberOfFiles << " files" << endl;
	cout << "Creating new file "<< outputFileName << endl;

	int nbOfEventsInEachFile = nbOfEvents / numberOfFiles;
	lcReader->close();
	lcReader->open( lcioFile );

	ostringstream ss;



	while( (evt = lcReader->readNextEvent()) != 0  ) {


		if ( nbOfWroteEvents == nbOfEventsInEachFile && nbOfWroteFiles != numberOfFiles ) {

			lcWriter->close();
			ss << nbOfWroteFiles;
			outputFileName = fileRoot + "_I"+ ss.str() +".slcio";
			lcWriter->open( outputFileName , LCIO::WRITE_NEW );

			cout << "Creating new file "<< outputFileName << endl;
			nbOfWroteFiles++;
			nbOfWroteEvents = 0;
			ss.str("");
		}

		lcWriter->writeEvent(evt);
		nbOfWroteEvents++;
	}


	lcReader->close();
	lcWriter->close();

	delete lcReader;
	delete lcWriter;

	return 0;
}

