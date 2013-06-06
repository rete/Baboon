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
	src/MarlinProcessor/BaboonProcessor.cc
	include/MarlinProcessor/BaboonProcessor.hh
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

# --------------------- Shower Splitter lib -------------------------------

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


# --------------------- Shower PDF lib -------------------------------

FILE (
	GLOB_RECURSE
	ShowerPDFProcessor_FILES
	include/MarlinProcessor/ShowerPDFProcessor.hh
	src/MarlinProcessor/ShowerPDFProcessor.cc
)
ADD_LIBRARY( ShowerPDFProcessor_lib SHARED ${ShowerPDFProcessor_FILES} )
SET_TARGET_PROPERTIES( ShowerPDFProcessor_lib PROPERTIES OUTPUT_NAME ShowerPDFProcessor )
INSTALL(
	TARGETS ShowerPDFProcessor_lib
	LIBRARY DESTINATION lib 
)

TARGET_LINK_LIBRARIES( 
	ShowerPDFProcessor_lib 
	${ROOT_COMPONENT_LIBRARIES} 
	Baboon 
)


FIND_LIBRARY( TSPECTRUM_LIB libSpectrum.so )



# --------------------- Overlay Estimator lib processor -------------------------------

FILE (
	GLOB_RECURSE
	OverlayEstimatorProcessor_FILES
	include/MarlinProcessor/OverlayEstimatorProcessor.hh
	src/MarlinProcessor/OverlayEstimatorProcessor.cc
)
ADD_LIBRARY( OverlayEstimator_lib SHARED ${OverlayEstimatorProcessor_FILES} )
SET_TARGET_PROPERTIES( OverlayEstimator_lib PROPERTIES OUTPUT_NAME OverlayEstimator )
INSTALL(
	TARGETS OverlayEstimator_lib
	LIBRARY DESTINATION lib 
)

TARGET_LINK_LIBRARIES( 
	OverlayEstimator_lib 
	${ROOT_COMPONENT_LIBRARIES} 
	Baboon
	${TSPECTRUM_LIB}
)



# --------------------- TestBeamCutsProcessor lib -------------------------------

FILE (
	GLOB_RECURSE
	TestBeamCutsProcessor_FILES
	include/MarlinProcessor/TestBeamCutsProcessor.hh
	src/MarlinProcessor/TestBeamCutsProcessor.cc
)
ADD_LIBRARY( TestBeamCutsProcessor_lib SHARED ${TestBeamCutsProcessor_FILES} )
SET_TARGET_PROPERTIES( TestBeamCutsProcessor_lib PROPERTIES OUTPUT_NAME TestBeamCutsProcessor )
INSTALL(
	TARGETS TestBeamCutsProcessor_lib
	LIBRARY DESTINATION lib 
)

TARGET_LINK_LIBRARIES( 
	TestBeamCutsProcessor_lib 
	${ROOT_COMPONENT_LIBRARIES} 
	Baboon 
)


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


# ---------------- Lcio split file ------------------
FILE( GLOB_RECURSE
	 LCIO_SPLIT_FILE_SRC_FILES 
	 main/LcioSplitFile.cc 
)
ADD_EXECUTABLE( LcioSplitFile ${LCIO_SPLIT_FILE_SRC_FILES} )
TARGET_LINK_LIBRARIES( 
	LcioSplitFile
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
#SET( 
#	OVERLAY_ESTIMATOR_SRC_FILES 
#	main/OverlayEstimator.cc 
#)
#ADD_EXECUTABLE( OverlayEstimator ${OVERLAY_ESTIMATOR_SRC_FILES} )
#TARGET_LINK_LIBRARIES( 
#	OverlayEstimator 
#	${ROOT_COMPONENT_LIBRARIES} 
#	Baboon
#	${TSPECTRUM_LIB}
#)


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


# ------------------ Overlay Reco Energy --------------------
SET( 
	OVERLAY_RECO_ENERGY_SRC_FILES 
	main/OverlayRecoEnergy.cc 
)
ADD_EXECUTABLE( OverlayRecoEnergy ${OVERLAY_RECO_ENERGY_SRC_FILES} )
TARGET_LINK_LIBRARIES( 
	OverlayRecoEnergy 
	${ROOT_COMPONENT_LIBRARIES} 
	Baboon
)


IF( BABOON_GUI )
	ADD_SUBDIRECTORY( ./gui )
ENDIF()



#-------------------- install ----------------------- 
INSTALL( 
	TARGETS 
	OverlayEvent LcioToRoot DrawEvent OverlayEstimatorAnalysis OverlayRecoEnergy LcioSplitFile
	DESTINATION 
	bin
)

