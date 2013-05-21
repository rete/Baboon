#################################################
# cmake file for building Baboon example processor
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
FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( LCIO REQUIRED )
FIND_PACKAGE( ROOT 5.32 REQUIRED )

FOREACH( pkg LCIO Marlin Baboon ROOT )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()

# include directories
INCLUDE_DIRECTORIES( include )

# require proper c++
ADD_DEFINITIONS( "-pedantic -Wunused-value -O2" )
ADD_DEFINITIONS("-Wno-long-long")

# set the processor files. Add other files for your binary in the list.
SET( 
	MY_SRC_FILES 
	myProc.cc 
	myProc.hh 
)

# build a shared library for marlin with src files
ADD_LIBRARY( 
	myProc 
	SHARED 
	${MY_SRC_FILES} 
)

# if your target is an executable, uncomment these lines and 
# adapt them to contain your 'main' and the needed src files
#ADD_EXECUTABLE(
#	myBin
#	${MY_SRC_FILES}
#)


# install the shared library in lib
INSTALL(
	TARGETS myProc
	LIBRARY DESTINATION lib 
)

# install the binary if there is one. Uncomment to use it. 
#INSTALL(
#	TARGETS myBin
#	RUNTIME DESTINATION bin 
#)

# display some variables and write them to cache
DISPLAY_STD_VARIABLES()
