#######################################################
# cmake file for building Baboon package
# @author rete , IPNL
#######################################################

# include directories
INCLUDE_DIRECTORIES( include )
INCLUDE_DIRECTORIES( BEFORE "${PROJECT_BINARY_DIR}" ) # generated BaboonConfig.h


# require proper c++
ADD_DEFINITIONS( "-pedantic -Wunused-value -O2 -std=c++11" )
ADD_DEFINITIONS("-Wno-long-long -Wreturn-type")
IF( OPENMP_FOUND )
	ADD_DEFINITIONS("-fopenmp")
ENDIF()



# create library
FOREACH( module
	Algorithm 
	Analysis
	Detector
	Geometry 
	Managers 
	Objects 
	Overlay 
	Reconstruction 
	Config 
	CfgParser  
	Utilities
	Monitoring
	Xml )
	
	FILE (
		GLOB_RECURSE
		BABOON_LIBRARY_FILES
		${BABOON_LIBRARY_FILES}
		include/${module}/*
		src/${module}/*
	)
ENDFOREACH()

ADD_SHARED_LIBRARY( Baboon ${BABOON_LIBRARY_FILES} )

TARGET_LINK_LIBRARIES( 
	Baboon
	${ROOT_COMPONENT_LIBRARIES} 
)

INSTALL(
	TARGETS Baboon
	LIBRARY DESTINATION lib 
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



IF( BABOON_GUI )
	#ADD_SUBDIRECTORY( ./gui )
	
	INCLUDE( ${ROOT_DICT_MACROS_FILE} )
	
	FILE (
		GLOB_RECURSE
		GUI_FILES
		include/Gui/*
		src/Gui/*
	)
	
	# include directories to be passed over to rootcint
	SET( ROOT_DICT_INCLUDE_DIRS
	    include
	    src
	)
	
	# get list of ALL Baboon GUI headers in include
	FILE( GLOB_RECURSE GUI_CXX_HEADERS_LIST include/Gui/EventNavigator.hh )
	
	FOREACH( mod 
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
		Xml )
		
		SET(
			ROOT_DICT_INPUT_HEADERS
			${ROOT_DICT_INPUT_HEADERS}
			./include/${mod}
		)
	ENDFOREACH()
	
	
	SET( 
		ROOT_DICT_INCLUDE_DIRS
		#${ROOT_DICT_INCLUDE_DIRS}
		${LCIO_INCLUDE_DIRS}
		${ROOT_INCLUDE_DIRS}
		${CLHEP_INCLUDE_DIRS}
		${Marlin_INCLUDE_DIRS}
		include/Gui
	    src/Gui
		./include
	)


	GEN_ROOT_DICT_LINKDEF_HEADER( baboon ${GUI_CXX_HEADERS_LIST} )
	GEN_ROOT_DICT_SOURCES( GUI.cxx )
	ADD_SHARED_LIBRARY( guiDict ${ROOT_DICT_OUTPUT_SOURCES} )
	INCLUDE_DIRECTORIES( "${ROOT_DICT_OUTPUT_DIR}" )
	
	ADD_CUSTOM_TARGET( dict DEPENDS guiDict )
	INSTALL_SHARED_LIBRARY( guiDict DESTINATION lib )
	
	ADD_SHARED_LIBRARY( BaboonGui SHARED ${GUI_FILES} )
	TARGET_LINK_LIBRARIES( guiDict Baboon BaboonGui )
	
	INSTALL(
		TARGETS BaboonGui
		LIBRARY DESTINATION lib 
	)
	SET( ALL_LIBS ${ALL_LIBS} ${mod} )
	TARGET_LINK_LIBRARIES( 
		BaboonGui 
		${ROOT_COMPONENT_LIBRARIES} 
	)
	
	
	
	# ------------------ Baboon GUI --------------------
	SET( 
		BABOON_GUI_SRC_FILES
		./src/MarlinProcessor/CaloHitCreator.cc
		./include/MarlinProcessor/CaloHitCreator.hh
		./main/Baboon.cc 
		${ROOT_DICT_OUTPUT_SOURCES}
	)
	ADD_EXECUTABLE( BaboonGUI ${BABOON_GUI_SRC_FILES} )
	SET_TARGET_PROPERTIES( BaboonGUI PROPERTIES OUTPUT_NAME Baboon )
	TARGET_LINK_LIBRARIES( 
		BaboonGUI 
		guiDict
		BaboonGui
		Baboon
	)
	
	
	#-------------------- install ----------------------- 
	INSTALL( 
		TARGETS 
		BaboonGUI
		DESTINATION 
		bin
	)
	
ENDIF()



#-------------------- install ----------------------- 
INSTALL( 
	TARGETS 
	OverlayEvent LcioToRoot LcioSplitFile
	DESTINATION 
	bin
)

