#######################################################
# cmake file for building Baboon package
# @author rete , IPNL
#######################################################

# include directories
INCLUDE_DIRECTORIES( include )
INCLUDE_DIRECTORIES( BEFORE "${PROJECT_BINARY_DIR}" ) # generated BaboonConfig.h


# require proper c++
ADD_DEFINITIONS( "-pedantic -Wunused-value -O2" )
ADD_DEFINITIONS("-Wno-long-long")
IF( OPENMP_FOUND )
	ADD_DEFINITIONS("-fopenmp")
ENDIF()

# create library
FOREACH( module
	Algorithm 
	Analysis 
	Geometry 
	Managers 
	Objects 
	Overlay 
	Reconstruction 
	Config 
	CfgParser 
	Utils 
	Utilities 
	Xml 
	GraphicalProcessor )
	
	FILE (
		GLOB_RECURSE
		BABOON_LIBRARY_FILES
		${BABOON_LIBRARY_FILES}
		include/${module}/*
		src/${module}/*
	)
ENDFOREACH()

SET( 
	BABOON_LIBRARY_FILES 
	${BABOON_LIBRARY_FILES} 
	src/MarlinProcessor/CaloHitCreator.cc
	include/MarlinProcessor/CaloHitCreator.hh
)

ADD_SHARED_LIBRARY( Baboon ${BABOON_LIBRARY_FILES} )
ADD_LIBRARY( Baboon_STATIC STATIC ${BABOON_LIBRARY_FILES} )

TARGET_LINK_LIBRARIES( 
	Baboon
	${ROOT_COMPONENT_LIBRARIES} 
)
TARGET_LINK_LIBRARIES( 
	Baboon_STATIC
	${ROOT_COMPONENT_LIBRARIES} 
)


SET_TARGET_PROPERTIES( Baboon_STATIC PROPERTIES OUTPUT_NAME Baboon )
INSTALL(
	TARGETS Baboon_STATIC Baboon
	ARCHIVE DESTINATION lib 
	LIBRARY DESTINATION lib 
)

FILE (
	GLOB_RECURSE
	ShowerSplitterProcessor_FILES
	include/MarlinProcessor/ShowerSplitterProcessor.hh
	src/MarlinProcessor/ShowerSplitterProcessor.cc
)
ADD_LIBRARY( ShowerSpliter_lib SHARED ${ShowerSplitterProcessor_FILES} )
SET_TARGET_PROPERTIES( ShowerSpliter_lib PROPERTIES OUTPUT_NAME ShowerSplitter )
INSTALL(
	TARGETS ShowerSpliter_lib
	LIBRARY DESTINATION lib 
)

TARGET_LINK_LIBRARIES( 
	ShowerSpliter_lib 
	${ROOT_COMPONENT_LIBRARIES} 
	Baboon 
)


FILE (
	GLOB_RECURSE
	GShowerSplitterProcessor_FILES
	include/MarlinProcessor/GShowerSplitterProcessor.hh
	src/MarlinProcessor/GShowerSplitterProcessor.cc
)
ADD_LIBRARY( GShowerSpliter_lib SHARED ${GShowerSplitterProcessor_FILES} )
SET_TARGET_PROPERTIES( GShowerSpliter_lib PROPERTIES OUTPUT_NAME GShowerSplitter )
INSTALL(
	TARGETS GShowerSpliter_lib
	LIBRARY DESTINATION lib 
)

TARGET_LINK_LIBRARIES( 
	GShowerSpliter_lib 
	${ROOT_COMPONENT_LIBRARIES} 
	Baboon 
)


FIND_LIBRARY( TSPECTRUM_LIB libSpectrum.so )


# ---------------- Overlay Event ------------------
FILE( GLOB_RECURSE
	 OVERLAY_EVENT_SRC_FILES 
	 main/OverlayEvent.cc 
)
ADD_EXECUTABLE( OverlayEvent ${OVERLAY_EVENT_SRC_FILES} )
TARGET_LINK_LIBRARIES( 
	OverlayEvent
	Baboon
)



# ------------------ LcioToRoot ---------------------
ADD_EXECUTABLE( LcioToRoot main/LcioToRoot.cc )
TARGET_LINK_LIBRARIES( LcioToRoot ${ROOT_COMPONENT_LIBRARIES} )



# ------------------ Draw event --------------------
SET( 
	DRAW_EVENT_SRC_FILES 
	main/DrawEvent.cc 
)
ADD_EXECUTABLE( DrawEvent ${DRAW_EVENT_SRC_FILES} )
TARGET_LINK_LIBRARIES( DrawEvent ${ROOT_COMPONENT_LIBRARIES} Baboon )


# ------------------ Overlay Estimator --------------------
SET( 
	OVERLAY_ESTIMATOR_SRC_FILES 
	main/OverlayEstimator.cc 
)
ADD_EXECUTABLE( OverlayEstimator ${OVERLAY_ESTIMATOR_SRC_FILES} )
TARGET_LINK_LIBRARIES( 
	OverlayEstimator 
	${ROOT_COMPONENT_LIBRARIES} 
	Baboon
	${TSPECTRUM_LIB}
)


# ------------------ Overlay Estimator Analysis --------------------
SET( 
	OVERLAY_ESTIMATOR_ANALYSIS_SRC_FILES 
	main/OverlayEstimatorAnalysis.cc 
)
ADD_EXECUTABLE( OverlayEstimatorAnalysis ${OVERLAY_ESTIMATOR_ANALYSIS_SRC_FILES} )
TARGET_LINK_LIBRARIES( 
	OverlayEstimatorAnalysis 
	${ROOT_COMPONENT_LIBRARIES} 
	Baboon
)



IF( BABOON_GUI )
	ADD_SUBDIRECTORY( ./gui )
ENDIF()



#-------------------- install ----------------------- 
INSTALL( 
	TARGETS 
	OverlayEvent LcioToRoot DrawEvent OverlayEstimator OverlayEstimatorAnalysis
	DESTINATION 
	bin
)

