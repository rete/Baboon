#################################################
# cmake file for building Baboon binaries
# @author Eté Rémi ,IPNL
# Copyright (c) CNRS / IPNL
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 FATAL_ERROR )
#################################################

# project name
PROJECT( Baboon )


SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 1 )


# --------------- ilcsoft dependencies -----------
LIST( APPEND CMAKE_MODULE_PATH $ILCSOFT )

# load default settings from ILCSOFT_CMAKE_MODULES
FIND_PACKAGE( ILCUTIL COMPONENTS ILCSOFT_CMAKE_MODULES REQUIRED )
INCLUDE( ilcsoft_default_settings )

# -------------- packages -------------
FIND_PACKAGE( LCIO REQUIRED )
FIND_PACKAGE( Marlin 1.0 REQUIRED )
FIND_PACKAGE( CLHEP REQUIRED )
FIND_PACKAGE( ROOT 5.27.06 REQUIRED COMPONENTS Gui Geom Eve Gdml TreePlayer RGL Ged EG GeomPainter GeomBuilder) # look for ROOT versions >= 5.27.06

FIND_PACKAGE( OpenMP)
if(OPENMP_FOUND)
	message("OPENMP FOUND")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

FOREACH( pkg LCIO ROOT Marlin CLHEP )
    IF( ${pkg}_FOUND )
        INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
        LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
        ADD_DEFINITIONS ( ${${pkg}_DEFINITIONS} )
    ENDIF()
ENDFOREACH()


OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" ON )

IF( INSTALL_DOC )

    FIND_PACKAGE( Doxygen )

    IF( DOXYGEN_FOUND )
        ADD_SUBDIRECTORY( ./doc )
    ELSE()
        MESSAGE( SEND_ERROR "Could not find doxygen required to build documentation" )
        MESSAGE( "Please install doxygen or set INSTALL_DOC to OFF" )
        MESSAGE( "" )
    ENDIF()

ENDIF()


# -------- generate BaboonConfig.h --------------------------------------------
CONFIGURE_FILE( "${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}Config.h.cmake.in"
    "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.h" @ONLY )



INSTALL( FILES "${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.h" DESTINATION "cmake" )


ADD_SUBDIRECTORY( ./source )


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()


# generate and install following configuration files
GENERATE_PACKAGE_CONFIGURATION_FILES( ${PROJECT_NAME}Config.cmake ${PROJECT_NAME}ConfigVersion.cmake ${PROJECT_NAME}LibDeps.cmake )

