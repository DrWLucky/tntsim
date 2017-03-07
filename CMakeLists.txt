# -*- mode: perl; -*-
# $Id: CMakeLists.txt 86065 2014-11-07 08:51:15Z gcosmo $

#--------------------CMake file for Geant 4.9.5.p01 with ROOT libs -----------
#Author: D. Perez-Loureiro 08/06/2012
#Modified: Brian Roeder 10/30/2013 for geant4.9.6, Margot_EU simulation
#Modified: Brian Roeder 10/08/2014 for geant4.10.00, Margot_EU simulation
#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(TNTSIM)

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
# Cmake Path for installation
#
#CMAKE_INSTALL_PREFIX 
#----------------------------------------------------------------------------
# Load some basic macros which are needed later on
EXECUTE_PROCESS(COMMAND root-config --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")

#---------------------------------------------------------------------------
# Create the directory for the ROOT files
#
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/root_files")

#
# Preprocessor derective for project source directory
#
add_definitions(-DTNTSIM_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")


#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
# Setup ROOT include directories
include(${Geant4_USE_FILE})

include_directories(${PROJECT_SOURCE_DIR}/include)

include_directories(${ROOT_INCLUDE_DIR})

#----------------------------------------------------------------------------
#Setup ROOT linking directories
#
link_directories( ${ROOT_LIBRARY_DIR})

#----------------------------------------------------------------------------
#Setup of compiler flags
if (CMAKE_SYSTEM_NAME MATCHES Darwin)
   EXEC_PROGRAM("sw_vers -productVersion | cut -d . -f 1-2" OUTPUT_VARIABLE MAC_OS_VERSION)
   MESSAGE("--- Found a Mac OS X System ${MAC_OS_VERSION}")
   if (CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE("--- Found GNU compiler collection")

      STRING(COMPARE EQUAL "10.5" "${MAC_OS_VERSION}" MAC_OS_10_5)
      IF(MAC_OS_10_5)
        SET(CMAKE_CXX_FLAGS "-m64")
        SET(CMAKE_Fortran_FLAGS "-m64")
      ENDIF(MAC_OS_10_5)

      SET(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")
      SET(CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_CXX_FLAGS} -flat_namespace -single_module -undefined dynamic_lookup")

      # Select flags.
      set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
      set(CMAKE_CXX_FLAGS_RELEASE        "-O2")
      set(CMAKE_CXX_FLAGS_DEBUG          "-g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_CXX_FLAGS_DEBUGFULL      "-g3 -fno-inline -Wnon-virtual-dtor -Wno-long-long -ansi -Wundef -Wcast-align -Wchar-subscripts -Wall -W -Wpointer-arith -Wformat-security -fno-exceptions -fno-check-new -fno-common")
      set(CMAKE_CXX_FLAGS_PROFILE        "-g3 -fno-inline -ftest-coverage -fprofile-arcs")
      set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g")
      set(CMAKE_C_FLAGS_RELEASE          "-O2")
      set(CMAKE_C_FLAGS_DEBUG            "-g -O2 -fno-reorder-blocks -fno-schedule-insns -fno-inline")
      set(CMAKE_C_FLAGS_DEBUGFULL        "-g3 -fno-inline -Wno-long-long -std=iso9899:1990 -Wundef -Wcast-align -Werror-implicit-function-declaration -Wchar-subscripts -Wall -W -Wpointer-arith -Wwrite-strings -Wformat-security -Wmissing-format-attribute -fno-common")
      set(CMAKE_C_FLAGS_PROFILE          "-g3 -fno-inline -ftest-coverage -fprofile-arcs")

      
      
      endif (CMAKE_COMPILER_IS_GNUCXX)

endif (CMAKE_SYSTEM_NAME MATCHES Darwin) 

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


#----------------------------------------------------------------------------
# Set everything needed for the root dictonary and create the
# dictionary
set(TNTSIM_HEADERS ${PROJECT_SOURCE_DIR}/include/TntActionInitialization.hh
                     ${PROJECT_SOURCE_DIR}/include/TntDetectorMessenger.hh
		     ${PROJECT_SOURCE_DIR}/include/TntEventMessenger.hh
		     ${PROJECT_SOURCE_DIR}/include/TntMuonPhysics.hh
		     ${PROJECT_SOURCE_DIR}/include/TntSteppingMessenger.hh
		     ${PROJECT_SOURCE_DIR}/include/TntTrajectory.hh
		     ${PROJECT_SOURCE_DIR}/include/TntWLSFiber.hh
		     ${PROJECT_SOURCE_DIR}/include/TntDataRecordTree.hh
		     ${PROJECT_SOURCE_DIR}/include/TntRecorderBase.hh
                     ${PROJECT_SOURCE_DIR}/include/TntDetectorMessenger.hh
		     ${PROJECT_SOURCE_DIR}/include/TntEMPhysics.hh
		     ${PROJECT_SOURCE_DIR}/include/TntGeneralPhysics.hh
		     ${PROJECT_SOURCE_DIR}/include/TntNuclearPhysics.hh
		     ${PROJECT_SOURCE_DIR}/include/TntPMTSD.hh
		     ${PROJECT_SOURCE_DIR}/include/TntRunAction.hh
		     ${PROJECT_SOURCE_DIR}/include/TntStackingAction.hh
                     ${PROJECT_SOURCE_DIR}/include/TntSteppingVerbose.hh
		     ${PROJECT_SOURCE_DIR}/include/TntUserEventInformation.hh
		     ${PROJECT_SOURCE_DIR}/include/TntWLSSlab.hh
		     ${PROJECT_SOURCE_DIR}/include/TntDetectorConstruction.hh
		     ${PROJECT_SOURCE_DIR}/include/TntEventAction.hh
		     ${PROJECT_SOURCE_DIR}/include/TntMainVolume.hh
		     ${PROJECT_SOURCE_DIR}/include/TntPhysicsList.hh
		     ${PROJECT_SOURCE_DIR}/include/TntPrimaryGeneratorAction.hh
		     ${PROJECT_SOURCE_DIR}/include/TntScintHit.hh
		     ${PROJECT_SOURCE_DIR}/include/TntScintSD.hh
		     ${PROJECT_SOURCE_DIR}/include/TntSteppingAction.hh
		     ${PROJECT_SOURCE_DIR}/include/TntTrackingAction.hh
		     ${PROJECT_SOURCE_DIR}/include/TntUserTrackInformation.hh
		     ${PROJECT_SOURCE_DIR}/include/menate_R.hh
		     )
#set(TNTSIM_LINKDEF ${PROJECT_SOURCE_DIR}/src/TNT2SIMbLinkDef.hh )
#set(TNTSIM_DICTIONARY TNT2SIMDict.cc) 
#ROOT_GENERATE_DICTIONARY("${TNT2SIM_HEADERS}" "${TNTSIM_LINKDEF}" "${TNTSIM_DICTIONARY}" "${ROOT_INCLUDE_DIR}")


#----------------------------------------------------------------------------
# Set the library version in the main CMakeLists.txt
SET(TNTSIM_MAJOR_VERSION 0)
SET(TNTSIM_MINOR_VERSION 0)
SET(TNTSIM_PATCH_VERSION 0)
SET(TNTSIM_VERSION
"${TNTSIM_MAJOR_VERSION}.${TNTSIM_MINOR_VERSION}.${TNTSIM_PATCH_VERSION}")
SET(TNTSIM_LIBRARY_PROPERTIES ${TNTSIM_LIBRARY_PROPERTIES}
    VERSION "${TNTSIM_VERSION}"
    SOVERSION "${TNTSIM_MAJOR_VERSION}"
    SUFFIX ".sl"
)





EXECUTE_PROCESS(COMMAND root-config --libs OUTPUT_VARIABLE ROOT_LD_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
set(CMAKE_EXE_LINKER_FLAGS ${ROOT_LD_FLAGS})
add_library(TNTSIM_lib SHARED tntsim.cc ${TNTSIM_DICTIONARY} ${sources} ${headers})
target_link_libraries(TNTSIM_lib ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
set_target_properties(TNTSIM_lib PROPERTIES ${TNTSIM_LIBRARY_PROPERTIES})

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(tntsim.exe tntsim.cc ${TNTSIM_DICTIONARY} ${sources} ${headers})
target_link_libraries(tntsim.exe ${Geant4_LIBRARIES} ${ROOT_LIBRARIES})
add_dependencies(tntsim.exe TNTSIM_lib)

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build TNTSIM. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
# BTR comment - Add here geant4 macros if used...
#
set(TNTSIM_SCRIPTS
  )

foreach(_script ${TNTSIM_SCRIPTS})
  configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()


#----------------------------------------------------------------------------
# For internal Geant4 use - but has no effect if you build this
# example standalone
#
add_custom_target(TNTSIM DEPENDS tntsim.exe)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
#install(TARGETS TNTSIM DESTINATION bin)

install(FILES ${headers}
        DESTINATION "./include"
        )
