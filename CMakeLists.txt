# This script generates CMakeLists.txt files for each example and adds them as external project
# The only thing necessary to build them against an existing OpenCMISS installation is to specify the 
# installation directory/prefix OPENCMISS_INSTALL_DIR containing the build context file (OpenCMISSBuildContext.cmake)
#
# You can change the build/install locations of the examples by providing
# OPENCMISS_EXAMPLES_BUILD_DIR and OPENCMISS_EXAMPLES_INSTALL_PREFIX.
#
# @author Daniel Wirtz @date April 9th, 2015
CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0-rc1 FATAL_ERROR)

# Paths and stuff
if (NOT OPENCMISS_INSTALL_DIR)
    set(OPENCMISS_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../install)
    get_filename_component(OPENCMISS_INSTALL_DIR ${OPENCMISS_INSTALL_DIR} ABSOLUTE)
endif()
if (NOT EXISTS ${OPENCMISS_INSTALL_DIR} OR NOT EXISTS ${OPENCMISS_INSTALL_DIR}/OpenCMISSBuildContext.cmake)
    message(FATAL_ERROR "OpenCMISS is not installed at '${OPENCMISS_INSTALL_DIR}'. Please specify OPENCMISS_INSTALL_DIR.")
endif()
include(${OPENCMISS_INSTALL_DIR}/OpenCMISSBuildContext.cmake)

project(OpenCMISS-Examples VERSION 1.0 LANGUAGES C)
include(ExternalProject)

# Set default build/install dirs unless specified
if (NOT OPENCMISS_EXAMPLES_BUILD_DIR)
    set(OPENCMISS_EXAMPLES_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
endif()
# Put the executable binaries into the opencmiss installation folder if not specified otherwise
#if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR NOT CMAKE_INSTALL_PREFIX)
#    set(CMAKE_INSTALL_PREFIX ${OPENCMISS_INSTALL_DIR}/examples CACHE STRING "Installation directory" FORCE)
#endif()
#message(STATUS "Installing examples to ${CMAKE_INSTALL_PREFIX}")

set(EXTRA_ARGS -DCMAKE_PREFIX_PATH=${OPENCMISS_PREFIX_PATH})
# For shared libs, use the correct install RPATH to enable binaries to find the shared libs.
# See http://www.cmake.org/Wiki/CMake_RPATH_handling
#if (BUILD_SHARED_LIBS) # always do it :-)
    list(APPEND EXTRA_ARGS 
        -DCMAKE_INSTALL_RPATH=${OPENCMISS_INSTALL_DIR}/lib
        -DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE
    )
#endif()

#get_filename_component(BINDIRNAME ${CMAKE_CURRENT_BINARY_DIR} NAME)
function(RECURSE curdir)
    file(GLOB content RELATIVE ${curdir} ${curdir}/*)
    SET(sources )
    foreach(entry ${content})
        set(filename ${curdir}/${entry})
        if(IS_DIRECTORY ${filename} AND (NOT filename STREQUAL CMAKE_CURRENT_BINARY_DIR))
                RECURSE(${filename})
        else()
            get_filename_component(fname ${filename} NAME)
            string(FIND "${curdir}" "${CMAKE_CURRENT_BINARY_DIR}" _TMP)
            if (fname STREQUAL Makefile 
                    AND (NOT curdir STREQUAL CMAKE_CURRENT_SOURCE_DIR)
                )
                
                # Add as external project right away
                GET_FULLNAME(${curdir} EXAMPLE_FULLNAME)
                message(STATUS "Configuring example ${EXAMPLE_FULLNAME}")
                
                #get_filename_component(EXAMPLE_NAME ${curdir} NAME)
                #string(TOLOWER EXAMPLE_NAME ${EXAMPLE_NAME})
                
                # Generation of CMakeLists.txt files from template
                #if (NOT EXISTS ${curdir}/CMakeLists.txt)
                    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.example.template.cmake
                        ${curdir}/CMakeLists.txt @ONLY)
                #endif()
                    
                ExternalProject_Add(${EXAMPLE_FULLNAME}
        		    PREFIX ${OPENCMISS_EXAMPLES_BUILD_DIR}
        		    TMP_DIR ${OPENCMISS_EXAMPLES_BUILD_DIR}/ep_tmp
        		    STAMP_DIR ${OPENCMISS_EXAMPLES_BUILD_DIR}/ep_stamps
        		
                    # Need empty download command, otherwise creation of external project fails with "no download info"
        		    DOWNLOAD_COMMAND ""
        		
        		    #--Configure step-------------
            		CMAKE_COMMAND ${CMAKE_COMMAND} -Wno-dev --no-warn-unused-cli # disables warnings for unused cmdline options
            		SOURCE_DIR ${curdir}
            		BINARY_DIR ${OPENCMISS_EXAMPLES_BUILD_DIR}/${EXAMPLE_FULLNAME}
            		# Need the prefix path in subprojects to find Iron and dependencies etc
                    # This variable is also set in the build context file.
            		CMAKE_ARGS 
            		    ${EXTRA_ARGS}
            		    -DCMAKE_INSTALL_PREFIX=${curdir}#${CMAKE_INSTALL_PREFIX}/${EXAMPLE_FULLNAME}
            	)
            endif()
        endif()
    endforeach()
endfunction()

# Create a suitable unique example target name from the example folder
function(GET_FULLNAME EXAMPLE_DIR VARNAME)
    string(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}" "" TMP "${EXAMPLE_DIR}")
    string(REGEX REPLACE "[/\\]" "_" TMP ${TMP})
    string(SUBSTRING ${TMP} 1 -1 TMP)
    #set(RES "")
    #while(NOT EXAMPLE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    #    get_filename_component(PARENT ${EXAMPLE_DIR} PATH)
    #    SET(EXAMPLE_DIR ${PARENT})
    #endwhile()
    set(${VARNAME} ${TMP} PARENT_SCOPE)
endfunction()

# Commence some real work here 
RECURSE(${CMAKE_CURRENT_SOURCE_DIR})