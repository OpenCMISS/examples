CMAKE_MINIMUM_REQUIRED(VERSION 3.2.0-rc1 FATAL_ERROR)
project(OpenCMISS-Examples VERSION 1.0 LANGUAGES C)
include(ExternalProject)

# Paths and stuff
if (NOT OPENCMISS_INSTALL_DIR)
    set(OPENCMISS_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../install)
    get_filename_component(OPENCMISS_INSTALL_DIR ${OPENCMISS_INSTALL_DIR} ABSOLUTE)
endif()
if (NOT EXISTS ${OPENCMISS_INSTALL_DIR} OR NOT EXISTS ${OPENCMISS_INSTALL_DIR}/iron-config.cmake)
    message(FATAL_ERROR "OpenCMISS is not installed at '${OPENCMISS_INSTALL_DIR}'. Please specify OPENCMISS_INSTALL_DIR.")
endif()
if (NOT OPENCMISS_EXAMPLES_BUILD_DIR)
    set(OPENCMISS_EXAMPLES_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR})
endif()
if (NOT OPENCMISS_EXAMPLES_INSTALL_PREFIX)
    set(OPENCMISS_EXAMPLES_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/../install)
endif()

list(APPEND EXAMPLE_EXTRA_DEFS -DCMAKE_PREFIX_PATH=${OPENCMISS_INSTALL_DIR})

function(RECURSE curdir)
    #message(STATUS "Entering ${curdir}")
    file(GLOB content RELATIVE ${curdir} ${curdir}/*)
    SET(sources )
    foreach(entry ${content})
        set(filename ${curdir}/${entry})
        if(IS_DIRECTORY ${filename})
            #list(FIND EXCLUDE_DIRS ${entry} EXCLUDE_ME)
            #if (EXCLUDE_ME EQUAL -1)
                RECURSE(${filename})
            #endif()
        else()
            get_filename_component(fname ${filename} NAME)
            # Process only c,h,cu files
            if (fname STREQUAL Makefile AND (NOT ${curdir} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR}))
                get_filename_component(EXAMPLE_NAME ${curdir} NAME)
                string(TOLOWER EXAMPLE_NAME ${EXAMPLE_NAME})
                # Generation of CMakeLists.txt files from template
                #if (NOT EXISTS ${curdir}/CMakeLists.txt)
                    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.example.template.cmake
                        ${curdir}/CMakeLists.txt @ONLY)
                #endif()
                    
                # Add as external project
                GET_FULLNAME(${curdir} EXAMPLE_FULLNAME)
                message(STATUS "Adding example ${EXAMPLE_FULLNAME}")
                ExternalProject_Add(${EXAMPLE_FULLNAME}
        		    PREFIX ${OPENCMISS_EXAMPLES_BUILD_DIR}
        		    TMP_DIR ${OPENCMISS_EXAMPLES_BUILD_DIR}/ep_tmp
        		    STAMP_DIR ${OPENCMISS_EXAMPLES_BUILD_DIR}/ep_stamps
        		
                    # Need empty download command, otherwise creation of external project fails with "no download info"
        		    DOWNLOAD_COMMAND ""
        		
        		    #--Configure step-------------
            		CMAKE_COMMAND ${CMAKE_COMMAND} --no-warn-unused-cli # disables warnings for unused cmdline options
            		SOURCE_DIR ${curdir}
            		BINARY_DIR ${OPENCMISS_EXAMPLES_BUILD_DIR}/${EXAMPLE_NAME}
            		CMAKE_ARGS ${EXAMPLE_EXTRA_DEFS}
            		
            		#--Build step-----------------
            		#SET(BUILD_CMD $)
                    #SET( INSTALL_CMD ${CMAKE_COMMAND} --build ${DIR} --target install)
            		#BUILD_COMMAND {CMAKE_COMMAND} --build ${DIR}
            		#--Install step---------------
            		# currently set as extra arg (above), somehow does not work
            		#INSTALL_DIR ${CMAKE_INSTALL_PREFIX} 
            		#INSTALL_COMMAND ${INSTALL_COMMAND}
            	)
            endif()
        endif()
    endforeach()
endfunction()

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
  
RECURSE(${CMAKE_CURRENT_SOURCE_DIR})