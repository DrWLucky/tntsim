# -*- mode: perl; -*-
# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(texansim)

## ---   Set the analysis package you want to use here --- ##
##
set(theAnalyzer "ROOT")
# set(theAnalyzer "XML")
# set(theAnalyzer "CSV")




#----------------------------------------------------------------------------
# Set some preprocessor definitions
#
    add_definitions(-DTEXAN_BUILD_DIR=\"${PROJECT_BINARY_DIR}\")
    add_definitions(-DTXS_MAX_HITS=10)

#----------------------------------------------------------------------------
# Define analysis package
#
if(theAnalyzer STREQUAL "ROOT")
    add_definitions(-DTEXAN_ANALYZE_ROOT)
elseif(theAnalyzer STREQUAL "XML")
    add_definitions(-DTEXAN_ANALYZE_XML)
elseif(theAnalyzer STREQUAL "CSV")
    add_definitions(-DTEXAN_ANALYZE_CSV)
endif(theAnalyzer STREQUAL "ROOT")


#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)
if(theAnalyzer STREQUAL "ROOT")
	EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
endif(theAnalyzer STREQUAL "ROOT")

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(texansim.exe texansim.cc ${sources} ${headers})
if(theAnalyzer STREQUAL "ROOT")
	EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
	set(CMAKE_EXE_LINKER_FLAGS ${ROOT_LD_FLAGS})
endif(theAnalyzer STREQUAL "ROOT")
target_link_libraries(texansim.exe ${Geant4_LIBRARIES})

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B1. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(TEXAN_SCRIPTS
    detectors.gdml
    )

foreach(_script ${TEXAN_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()

# add a target to generate API documentation with Doxygen

find_package(Doxygen)
    if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
    add_custom_target(doc
		      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		      COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
    endif(DOXYGEN_FOUND)


#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(texansim DEPENDS texansim.exe)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS texansim.exe DESTINATION bin)
