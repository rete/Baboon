#################################################
# cmake file for building Baboon binaries
# @author Eté Rémi ,IPNL
# Copyright (c) CNRS / IPNL
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
#################################################

# project name
PROJECT( myBaboonProc )


# --------------- ilcsoft dependencies -----------
LIST( APPEND CMAKE_MODULE_PATH $ILCSOFT )

# load default settings from ILCSOFT_CMAKE_MODULES
FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )
INCLUDE( ilcsoft_default_settings )

# -------------- packages -------------
FIND_PACKAGE( Baboon REQUIRED )

IF( Baboon_FOUND )
    INCLUDE_DIRECTORIES( ${Baboon_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${Baboon_LIBRARIES} )
    ADD_DEFINITIONS ( ${Baboon_DEFINITIONS} )
ENDIF()


INCLUDE_DIRECTORIES( include )


# require proper c++
ADD_DEFINITIONS( "-pedantic -Wunused-value -O2" )
ADD_DEFINITIONS("-Wno-long-long")

ADD_LIBRARY( 
	myProc 
	SHARED 
	src/myProc.cc 
	include/myProc.hh 
)

INSTALL(
	TARGETS myProc
	LIBRARY DESTINATION lib 
)


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
