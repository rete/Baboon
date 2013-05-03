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
FOREACH( mod Algorithm Analysis Geometry Managers Objects Overlay Reconstruction Config CfgParser Utils Utilities )
	FILE (
		GLOB_RECURSE
		${mod}_FILES
		include/${mod}/*
		src/${mod}/*
	)
	ADD_LIBRARY( ${mod} SHARED ${${mod}_FILES} )
	ADD_LIBRARY( ${mod}_STATIC STATIC ${${mod}_FILES} )
	SET_TARGET_PROPERTIES( ${mod}_STATIC PROPERTIES OUTPUT_NAME ${mod} )
	INSTALL(
		TARGETS ${mod}_STATIC ${mod}
		ARCHIVE DESTINATION lib 
		LIBRARY DESTINATION lib 
	)
ENDFOREACH()


TARGET_LINK_LIBRARIES( Utils Objects )
TARGET_LINK_LIBRARIES( Managers Analysis CfgParser )
TARGET_LINK_LIBRARIES( Config CfgParser )
TARGET_LINK_LIBRARIES( Algorithm Objects )


FIND_LIBRARY( TSPECTRUM_LIB libSpectrum.so )


# ---------------- Overlay Event ------------------
FILE( GLOB_RECURSE
	 OVERLAY_EVENT_SRC_FILES 
	 main/OverlayEvent.cc 
)
ADD_EXECUTABLE( OverlayEvent ${OVERLAY_EVENT_SRC_FILES} )
TARGET_LINK_LIBRARIES( 
	OverlayEvent  
	Managers
	Objects
	Overlay
	Utils
	Config
	CfgParser
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
TARGET_LINK_LIBRARIES( DrawEvent ${ROOT_COMPONENT_LIBRARIES} CfgParser Config )


# ------------------ Overlay Estimator --------------------
SET( 
	OVERLAY_ESTIMATOR_SRC_FILES 
	main/OverlayEstimator.cc 
)
ADD_EXECUTABLE( OverlayEstimator ${OVERLAY_ESTIMATOR_SRC_FILES} )
TARGET_LINK_LIBRARIES( 
	OverlayEstimator 
	${ROOT_COMPONENT_LIBRARIES} 
	CfgParser 
	Config 
	Analysis
	Objects
	Geometry
	Algorithm
	Managers
	${TSPECTRUM_LIB}
	Utils
	Utilities
)


# ------------------ Shower Spliter --------------------
SET( 
	SHOWER_SPLITER_SRC_FILES 
	main/ShowerSpliter.cc 
)
ADD_EXECUTABLE( ShowerSpliter ${SHOWER_SPLITER_SRC_FILES} )
TARGET_LINK_LIBRARIES( 
	ShowerSpliter 
	${ROOT_COMPONENT_LIBRARIES} 
	CfgParser 
	Config 
	Analysis
	Objects
	Geometry
	Algorithm
	Managers
	${TSPECTRUM_LIB}
	Utils
	Utilities
)



#-------------------- install ----------------------- 
INSTALL( 
	TARGETS 
	OverlayEvent LcioToRoot DrawEvent ShowerSpliter OverlayEstimator
	DESTINATION 
	bin
)
