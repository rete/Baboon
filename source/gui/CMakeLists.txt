#######################################################
# cmake file for building Baboon gui
# @author rete , IPNL
#######################################################

# include directories
INCLUDE_DIRECTORIES( include )
INCLUDE_DIRECTORIES( BEFORE "${PROJECT_BINARY_DIR}" ) # generated BaboonConfig.h


# require proper c++
ADD_DEFINITIONS( "-pedantic -Wunused-value -O2" )
ADD_DEFINITIONS("-Wno-long-long")

FILE (
	GLOB_RECURSE
	GUI_FILES
	include/*
	src/*
)

INCLUDE( ${ROOT_DICT_MACROS_FILE} )

# include directories to be passed over to rootcint
SET( ROOT_DICT_INCLUDE_DIRS
    include
    src
)

# get list of ALL Baboon GUI headers in include
FILE( GLOB_RECURSE GUI_CXX_HEADERS_LIST ../include/* include/EventNavigator.hh )


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
		../include/${mod}
	)
ENDFOREACH()
SET( 
	ROOT_DICT_INCLUDE_DIRS
	#${ROOT_DICT_INCLUDE_DIRS}
	${LCIO_INCLUDE_DIRS}
	${ROOT_INCLUDE_DIRS}
	${CLHEP_INCLUDE_DIRS}
	${Marlin_INCLUDE_DIRS}
	include
    src
	../include
)
MESSAGE( "ROOT_DICT_INCLUDE_DIRS : ${ROOT_DICT_INCLUDE_DIRS}")
MESSAGE( "includes : " ${GUI_CXX_HEADERS_LIST} )

GEN_ROOT_DICT_LINKDEF_HEADER( baboon ${GUI_CXX_HEADERS_LIST} )
GEN_ROOT_DICT_SOURCES( GUI.cxx )
#MESSAGE( "dict sources : ${ROOT_DICT_OUTPUT_SOURCES}" )
ADD_SHARED_LIBRARY( guiDict ${ROOT_DICT_OUTPUT_SOURCES} )

INCLUDE_DIRECTORIES( "${ROOT_DICT_OUTPUT_DIR}" )


ADD_SHARED_LIBRARY( BaboonGui SHARED ${GUI_FILES} )
TARGET_LINK_LIBRARIES( guiDict Baboon BaboonGui )

ADD_CUSTOM_TARGET( dict DEPENDS guiDict )
INSTALL_SHARED_LIBRARY( guiDict DESTINATION lib )

ADD_LIBRARY( BaboonGui_STATIC STATIC ${GUI_FILES} )
SET_TARGET_PROPERTIES( BaboonGui_STATIC PROPERTIES OUTPUT_NAME BaboonGui )
INSTALL(
	TARGETS BaboonGui_STATIC BaboonGui
	ARCHIVE DESTINATION lib 
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
	../src/MarlinProcessor/CaloHitCreator.cc
	../include/MarlinProcessor/CaloHitCreator.hh
	./BaboonGui.cc 
)
ADD_EXECUTABLE( BaboonGUI ${BABOON_GUI_SRC_FILES} )
TARGET_LINK_LIBRARIES( 
	BaboonGUI 
	#guiDict
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

