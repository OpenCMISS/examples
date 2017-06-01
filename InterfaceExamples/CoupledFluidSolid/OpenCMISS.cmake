###########################
# *DO NOT CHANGE THIS FILE*
###########################
#
# Prepares the use of OpenCMISS by scanning at various path locations for OpenCMISS (SDK) installations.
#
# Search order:
#     - OPENCMISS_INSTALL_DIR variable
#     - OPENCMISS_INSTALL_DIR environment variable
#     - OPENCMISS_SDK_DIR variable
#     - OPENCMISS_SDK_DIR environment variable
#     - [Windows only] The Registry is searched additionally for installed SDKs.
#
# If found
#     - Includes the toolchain config script of the found opencmiss installation to enable toolchain selection logic (non-windows)

# Convenience: The OPENCMISS_INSTALL_DIR may also be defined in the environment and various other places
if (DEFINED OPENCMISS_INSTALL_DIR)
    get_filename_component(OPENCMISS_INSTALL_DIR "${OPENCMISS_INSTALL_DIR}" ABSOLUTE)
    if (EXISTS "${OPENCMISS_INSTALL_DIR}")
        message(STATUS "Using specified installation directory '${OPENCMISS_INSTALL_DIR}'")
    else()
        message(WARNING "The specified OPENCMISS_INSTALL_DIR '${OPENCMISS_INSTALL_DIR}' does not exist. Skipping.")
        unset(OPENCMISS_INSTALL_DIR)
    endif()
endif()
if(NOT OPENCMISS_INSTALL_DIR AND NOT "$ENV{OPENCMISS_INSTALL_DIR}" STREQUAL "")
    file(TO_CMAKE_PATH "$ENV{OPENCMISS_INSTALL_DIR}" OPENCMISS_INSTALL_DIR)
    get_filename_component(OPENCMISS_INSTALL_DIR "${OPENCMISS_INSTALL_DIR}" ABSOLUTE)
    if (EXISTS "${OPENCMISS_INSTALL_DIR}")
        message(STATUS "Using environment installation directory '${OPENCMISS_INSTALL_DIR}'")
    else()
        message(WARNING "The environment variable OPENCMISS_INSTALL_DIR='${OPENCMISS_INSTALL_DIR}' contains an invalid path. Skipping.")
        unset(OPENCMISS_INSTALL_DIR)
    endif()
endif()
if (NOT OPENCMISS_INSTALL_DIR AND DEFINED OPENCMISS_SDK_DIR)
    get_filename_component(OPENCMISS_SDK_DIR "${OPENCMISS_SDK_DIR}" ABSOLUTE)
    if (EXISTS "${OPENCMISS_SDK_DIR}")
        message(STATUS "Using SDK installation directory: ${OPENCMISS_SDK_DIR}")
        set(OPENCMISS_INSTALL_DIR "${OPENCMISS_SDK_DIR}")
    else()
        message(WARNING "The specified OPENCMISS_SDK_DIR '${OPENCMISS_SDK_DIR}' does not exist. Skipping.")
    endif()
endif()    
if(NOT OPENCMISS_INSTALL_DIR AND NOT "$ENV{OPENCMISS_SDK_DIR}" STREQUAL "")
    file(TO_CMAKE_PATH "$ENV{OPENCMISS_SDK_DIR}" OPENCMISS_SDK_DIR)
    get_filename_component(OPENCMISS_SDK_DIR "${OPENCMISS_SDK_DIR}" ABSOLUTE)
    if (EXISTS "${OPENCMISS_SDK_DIR}")
        message(STATUS "Using environment SDK installation directory: ${OPENCMISS_SDK_DIR}")
        set(OPENCMISS_INSTALL_DIR "${OPENCMISS_SDK_DIR}")
    else()
        message(WARNING "The environment variable OPENCMISS_SDK_DIR='${OPENCMISS_SDK_DIR}' contains an invalid path. Skipping.")
        unset(OPENCMISS_SDK_DIR)
    endif()
endif()
# On windows: check the registry for installed OpenCMISS SDKs
if(NOT OPENCMISS_INSTALL_DIR AND WIN32)
    foreach(ROOT HKEY_CURRENT_USER HKEY_LOCAL_MACHINE)
        foreach(PACKAGE OpenCMISSUserSDK OpenCMISSDeveloperSDK)
            set(_REG_KEY "[${ROOT}\\Software\\Auckland Bioengineering Institute\\${PACKAGE}]")
            get_filename_component(OPENCMISS_INSTALL_DIR "${_REG_KEY}" ABSOLUTE)
            #message(STATUS "Trying registry key ${_REG_KEY}: ${OPENCMISS_INSTALL_DIR}")
            if (EXISTS "${OPENCMISS_INSTALL_DIR}")
                message(STATUS "Found ${PACKAGE} in Windows registry key ${_REG_KEY}: ${OPENCMISS_INSTALL_DIR}")
                break()
            else()
                unset(OPENCMISS_INSTALL_DIR)
            endif()         
        endforeach()
        if (EXISTS "${OPENCMISS_INSTALL_DIR}")
            break()
        endif()
    endforeach()
endif()

if(OPENCMISS_INSTALL_DIR)
    message(STATUS "OPENCMISS_INSTALL_DIR=${OPENCMISS_INSTALL_DIR}")
    # Use the OpenCMISS scripts to also allow choosing a separate toolchain
    # This file is located at the opencmiss installation rather than the local example
    # as it avoids file replication and makes maintenance much easier
    if (TOOLCHAIN)
        set(_OCTC ${OPENCMISS_INSTALL_DIR}/cmake/OCToolchainCompilers.cmake)
        if (EXISTS "${_OCTC}")
            include(${_OCTC})
        else()
            message(WARNING "TOOLCHAIN specified but OpenCMISS config script could not be found at ${_OCTC}. Trying CMake defaults.")
        endif()
        unset(_OCTC)
    endif()
    set(OpenCMISS_DIR "${OPENCMISS_INSTALL_DIR}" CACHE PATH "Path to the found OpenCMISS (SDK) installation")    
else()
    message(WARNING "No OpenCMISS (SDK) installation directory detected. Using default system environment.")
endif()
 
