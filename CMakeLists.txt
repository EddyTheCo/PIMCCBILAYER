# $Id: CMakeLists.txt 68058 2013-03-13 14:47:43Z gcosmo $
#
#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(PIMCC)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if(USEROOT)
    message( STATUS "This build is USING ROOT")
add_definitions(-DUSEROOT)
set(ROOT_INSTALL /home/camelud/GEANT4/ROOT/build)
include_directories(${ROOT_INSTALL}/include)
endif()


if(SAVECONF)
    message( STATUS "This build will Save CONF ")
add_definitions(-DSAVECONF)
endif()
if(UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -g -std=c++11 -O3 -march=native")        ## Optimize
    set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
endif()




include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_SOURCE_DIR}/src)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)


add_executable(PIMCC PIMCC.cc ${sources} ${headers}  )

if(USEROOT)
target_link_libraries(PIMCC -lstdc++ -L${ROOT_INSTALL}/lib  -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lm -ldl -rdynamic)
endif()
#if(CMAKE_COMPILER_IS_GNUCXX)
 #   set(CMAKE_CXX_FLAGS "-O3")        ## Optimize
  #  set(CMAKE_EXE_LINKER_FLAGS "-s")  ## Strip binary
#endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build B4a. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(SCRIPTS
#  INIMaterial
#  INISource
  input
  rootScript.cpp
 )

foreach(_script ${SCRIPTS})
 configure_file(
    ${PROJECT_SOURCE_DIR}/${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
file(COPY ${CMAKE_SOURCE_DIR}/TESTS DESTINATION ${PROJECT_BINARY_DIR})

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS PIMCC DESTINATION bin)
#else()
#   message("#ERROR need to install root from https://root.cern.ch ")
#endif()
