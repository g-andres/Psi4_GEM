# LibxcConfig.cmake
# ------------------
#
# Libxc cmake module.
# This module sets the following variables in your project:
#
# ::
#
#   Libxc_FOUND - true if Libxc and all required components found on the system
#   Libxc_VERSION - Libxc version in format Major.Minor.Release
#   Libxc_INCLUDE_DIRS - Directory where Libxc header is located.
#   Libxc_INCLUDE_DIR - same as DIRS
#   Libxc_DEFINITIONS: Definitions necessary to use Libxc, namely USING_Libxc.
#   Libxc_LIBRARIES - Libxc library to link against.
#   Libxc_LIBRARY - same as LIBRARIES
#
#
# Available components: shared static
#
# ::
#
#   shared - search for only shared library
#   static - search for only static library
#   C - search for C library only, even if Fortran available
#   Fortran - search for Fortran library (C library always present)
#
#
# Exported targets:
#
# ::
#
# If Libxc is found, this module defines the following :prop_tgt:`IMPORTED`
# targets. Target is shared _or_ static, so, for both, use separate, not
# overlapping, installations. ::
#
#   Libxc::xc - the main Libxc library with header & defs attached.
#   Libxc::xcf03 - the Fortran 2003 Libxc library that depends on Libxc::xc
#
#
# Suggested usage:
#
# ::
#
#   find_package(Libxc)
#   find_package(Libxc 5.1.0 EXACT CONFIG REQUIRED COMPONENTS shared C)
#
#
# The following variables can be set to guide the search for this package:
#
# ::
#
#   Libxc_DIR - CMake variable, set to directory containing this Config file
#   CMAKE_PREFIX_PATH - CMake variable, set to root directory of this package
#   PATH - environment variable, set to bin directory of this package
#   CMAKE_DISABLE_FIND_PACKAGE_Libxc - CMake variable, disables
#     find_package(Libxc) when not REQUIRED, perhaps to force internal build


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was LibxcConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

set(PN Libxc)
set (_valid_components
    static
    shared
    C
    Fortran
)

# check library style component
if(OFF)
    set(${PN}_shared_FOUND 1)
else()
    set(${PN}_static_FOUND 1)
endif()
list(FIND ${PN}_FIND_COMPONENTS "shared" _seek_shared)
list(FIND ${PN}_FIND_COMPONENTS "static" _seek_static)

# check library language component
set(${PN}_C_FOUND 1)
list(FIND ${PN}_FIND_COMPONENTS "C" _seek_C)
if(OFF)
    set(${PN}_Fortran_FOUND 1)
endif()
list(FIND ${PN}_FIND_COMPONENTS "Fortran" _seek_Fortran)

set(${PN}_DEFINITIONS USING_${PN})

check_required_components(${PN})

#-----------------------------------------------------------------------------
# Don't include targets if this file is being picked up by another
# project which has already built this as a subproject
#-----------------------------------------------------------------------------
if(NOT TARGET ${PN}::xc)

    if(_seek_Fortran GREATER -1)
        set(_target "xcf03")
        include("${CMAKE_CURRENT_LIST_DIR}/${PN}Targets-Fortran.cmake")
    elseif(_seek_C GREATER -1)
        set(_target "xc")
        include("${CMAKE_CURRENT_LIST_DIR}/${PN}Targets-C.cmake")
    elseif(OFF)
        set(_target "xcf03")
        include("${CMAKE_CURRENT_LIST_DIR}/${PN}Targets-Fortran.cmake")
    else()
        set(_target "xc")
        include("${CMAKE_CURRENT_LIST_DIR}/${PN}Targets-C.cmake")
    endif()

    get_property(_loc TARGET ${PN}::${_target} PROPERTY LOCATION)
    set(${PN}_LIBRARY ${_loc})
    get_property(_ill TARGET ${PN}::${_target} PROPERTY INTERFACE_LINK_LIBRARIES)
    set(${PN}_LIBRARIES ${_ill})

    get_property(_id TARGET ${PN}::${_target} PROPERTY INCLUDE_DIRECTORIES)
    set(${PN}_INCLUDE_DIR ${_id})
    get_property(_iid TARGET ${PN}::${_target} PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    set(${PN}_INCLUDE_DIRS ${_iid})

    #message("Libxc::${_target}")
    #message("loc ${_loc}")
    #message("ill ${_ill}")
    #message("id  ${_id}")
    #message("iid ${_iid}")
endif()
