# Libint2Config.cmake
# -------------------
#
# Libint2 cmake module.
# This module sets the following variables in your project:
#
# ::
#
#   Libint2_FOUND - true if Libint2 and all required components found on the system
#   Libint2_VERSION - Libint2 version in format Major.Minor.Release
#   Libint2_INCLUDE_DIRS - Directories where Libint2 and libderiv headers are located.
#   Libint2_INCLUDE_DIR - same as DIRS
#   Libint2_LIBRARIES - Libint2 and libderiv libraries to link against.
#   Libint2_LIBRARY - same as LIBRARIES
#   Libint2_MAX_AM_ERI - maximum angular momentum level of Libint2 libraries
#
#
# Available components:
#
# ::
#
#   shared - search for only shared library
#   static - search for only static library
#
#   e[2, 8] - search for library including energy ERI integrals with angular momentum >= this integer
#   g[2, 8] - search for library including gradient ERI integrals with angular momentum >= this integer
#   h[2, 8] - search for library including Hessian ERI integrals with angular momentum >= this integer
#   eri2_e[2, 8] - search for library including energy ERI2 integrals with angular momentum >= this integer
#   eri2_g[2, 8] - search for library including gradient ERI2 integrals with angular momentum >= this integer
#   eri2_h[2, 8] - search for library including Hessian ERI2 integrals with angular momentum >= this integer
#   eri3_e[2, 8] - search for library including energy ERI2 integrals with angular momentum >= this integer
#   eri3_g[2, 8] - search for library including gradient ERI2 integrals with angular momentum >= this integer
#   eri3_h[2, 8] - search for library including Hessian ERI2 integrals with angular momentum >= this integer
#
#                    sph        cart       shell_set  used_by
#                    --------   --------   ---------  -------
#   sss - search for standard + standard + standard = mpqc4
#   sso - search for                     + orca
#   sis - search for          + intv3    + standard
#   sio - search for                     + orca
#   sgs - search for          + gamess   + standard = gamess
#   sgo - search for                     + orca
#   sos - search for          + orca     + standard
#   soo - search for                     + orca     = orca
#   sbs - search for          + bagel    + standard = bagel
#   sbo - search for                     + orca
#   gss - search for gaussian + standard + standard = psi4
#   gso - search for                     + orca
#   gis - search for          + intv3    + standard
#   gio - search for                     + orca
#   ggs - search for          + gamess   + standard
#   ggo - search for                     + orca
#   gos - search for          + orca     + standard
#   goo - search for                     + orca
#   gbs - search for          + bagel    + standard
#   gbo - search for                     + orca
#
#
# Exported targets:
#
# ::
#
# If Libint2 is found, this module defines the following :prop_tgt:`IMPORTED`
# targets. ::
#
##   Libint2::int2 - library only
#   Libint2::cxx - library + C++11 API
#
#
# Suggested usage:
#
# ::
#
#   find_package(Libint2)
#   find_package(Libint2 2.7.0 CONFIG REQUIRED COMPONENTS shared gss e5 g5)
#
#
# The following variables can be set to guide the search for this package:
#
# ::
#
#   Libint2_DIR - CMake variable, set to directory containing this Config file
#   CMAKE_PREFIX_PATH - CMake variable, set to root directory of this package
#   PATH - environment variable, set to bin directory of this package
#   CMAKE_DISABLE_FIND_PACKAGE_Libint2 - CMake variable, disables
#     find_package(Libint2) when not REQUIRED, perhaps to force internal build


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was Libint2Config.cmake.in                            ########

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

set(PN Libint2)
set (_valid_components
    static
    shared

    e2
    e3
    e4
    e5
    e6
    e7
    e8
    g2
    g3
    g4
    g5
    g6
    g7
    g8
    h2
    h3
    h4
    h5
    h6
    h7
    h8

    eri2_e2
    eri2_e3
    eri2_e4
    eri2_e5
    eri2_e6
    eri2_e7
    eri2_e8
    eri2_g2
    eri2_g3
    eri2_g4
    eri2_g5
    eri2_g6
    eri2_g7
    eri2_g8
    eri2_h2
    eri2_h3
    eri2_h4
    eri2_h5
    eri2_h6
    eri2_h7
    eri2_h8

    eri3_e2
    eri3_e3
    eri3_e4
    eri3_e5
    eri3_e6
    eri3_e7
    eri3_e8
    eri3_g2
    eri3_g3
    eri3_g4
    eri3_g5
    eri3_g6
    eri3_g7
    eri3_g8
    eri3_h2
    eri3_h3
    eri3_h4
    eri3_h5
    eri3_h6
    eri3_h7
    eri3_h8

    sss
    sso
    sis
    sio
    sgs
    sgo
    sos
    soo
    sbs
    sbo
    gss
    gso
    gis
    gio
    ggs
    ggo
    gos
    goo
    gbs
    gbo
)

# check library style component
if(OFF)
    set(${PN}_shared_FOUND 1)
endif()
if(ON)
    set(${PN}_static_FOUND 1)
endif()
list(FIND ${PN}_FIND_COMPONENTS "shared" _seek_shared)
list(FIND ${PN}_FIND_COMPONENTS "static" _seek_static)

# check AM & derivative component
set(${PN}_MAX_AM_ERI 5)
foreach(_eri e2;e3;e4;e5;g2;g3;g4;h2;h3;eri2_e2;eri2_e3;eri2_e4;eri2_e5;eri2_e6;eri2_g2;eri2_g3;eri2_g4;eri2_g5;eri2_h2;eri2_h3;eri2_h4;eri3_e2;eri3_e3;eri3_e4;eri3_e5;eri3_e6;eri3_g2;eri3_g3;eri3_g4;eri3_g5;eri3_h2;eri3_h3;eri3_h4)
    set(${PN}_${_eri}_FOUND 1)
endforeach()

# check orderings component
if    ((2 EQUAL 1) AND (1 EQUAL 1) AND (1 EQUAL 1))
    set(_ordering "sss")
elseif((2 EQUAL 1) AND (1 EQUAL 1) AND (1 EQUAL 2))
    set(_ordering "sso")
elseif((2 EQUAL 1) AND (1 EQUAL 2) AND (1 EQUAL 1))
    set(_ordering "sis")
elseif((2 EQUAL 1) AND (1 EQUAL 2) AND (1 EQUAL 2))
    set(_ordering "sio")
elseif((2 EQUAL 1) AND (1 EQUAL 3) AND (1 EQUAL 1))
    set(_ordering "sgs")
elseif((2 EQUAL 1) AND (1 EQUAL 3) AND (1 EQUAL 2))
    set(_ordering "sgo")
elseif((2 EQUAL 1) AND (1 EQUAL 4) AND (1 EQUAL 1))
    set(_ordering "sos")
elseif((2 EQUAL 1) AND (1 EQUAL 4) AND (1 EQUAL 2))
    set(_ordering "soo")
elseif((2 EQUAL 1) AND (1 EQUAL 5) AND (1 EQUAL 1))
    set(_ordering "sbs")
elseif((2 EQUAL 1) AND (1 EQUAL 5) AND (1 EQUAL 2))
    set(_ordering "sbo")
elseif((2 EQUAL 2) AND (1 EQUAL 1) AND (1 EQUAL 1))
    set(_ordering "gss")
elseif((2 EQUAL 2) AND (1 EQUAL 1) AND (1 EQUAL 2))
    set(_ordering "gso")
elseif((2 EQUAL 2) AND (1 EQUAL 2) AND (1 EQUAL 1))
    set(_ordering "gis")
elseif((2 EQUAL 2) AND (1 EQUAL 2) AND (1 EQUAL 2))
    set(_ordering "gio")
elseif((2 EQUAL 2) AND (1 EQUAL 3) AND (1 EQUAL 1))
    set(_ordering "ggs")
elseif((2 EQUAL 2) AND (1 EQUAL 3) AND (1 EQUAL 2))
    set(_ordering "ggo")
elseif((2 EQUAL 2) AND (1 EQUAL 4) AND (1 EQUAL 1))
    set(_ordering "gos")
elseif((2 EQUAL 2) AND (1 EQUAL 4) AND (1 EQUAL 2))
    set(_ordering "goo")
elseif((2 EQUAL 2) AND (1 EQUAL 5) AND (1 EQUAL 1))
    set(_ordering "gbs")
elseif((2 EQUAL 2) AND (1 EQUAL 5) AND (1 EQUAL 2))
    set(_ordering "gbo")
else()
    message(STATUS "${PN}Config: indeterminate orderings")
endif()
set(${PN}_${_ordering}_FOUND 1)


# thanks, https://stackoverflow.com/a/9328525
function(dump_cmake_variables)
    get_cmake_property(_variableNames VARIABLES)
    list (SORT _variableNames)
    set(founds "")
    foreach (_variableName ${_variableNames})
        if (ARGV0)
            unset(MATCHED)
            string(REGEX MATCH ${ARGV0} MATCHED ${_variableName})
            if (NOT MATCHED)
                continue()
            endif()
            if (NOT ${${_variableName}})
                continue()
            endif()
        endif()
        list(APPEND found ${CMAKE_MATCH_1})
    endforeach()
    message(STATUS "${ARGV1}${found}")
endfunction()


if(NOT CMAKE_REQUIRED_QUIET)
    message(STATUS "${PN}Config components requested: ${${PN}_FIND_COMPONENTS}")
    dump_cmake_variables("^Libint2_([A-Za-z0-9_]+)_FOUND$" "${PN}Config components found: ")
endif()

check_required_components(${PN})

# make detectable the various cmake modules exported alongside
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})

#-----------------------------------------------------------------------------
# Don't include targets if this file is being picked up by another
# project which has already built this as a subproject
#-----------------------------------------------------------------------------
if(NOT TARGET ${PN}::cxx)
    if(_seek_static GREATER -1)
        include("${CMAKE_CURRENT_LIST_DIR}/${PN}Targets-static.cmake")
    elseif(_seek_shared GREATER -1)
        include("${CMAKE_CURRENT_LIST_DIR}/${PN}Targets-shared.cmake")
    elseif(OFF)
        include("${CMAKE_CURRENT_LIST_DIR}/${PN}Targets-shared.cmake")
    elseif(ON)
        include("${CMAKE_CURRENT_LIST_DIR}/${PN}Targets-static.cmake")
    endif()

    include(CMakeFindDependencyMacro)
    if(NOT TARGET Eigen3::Eigen)
        find_dependency(Eigen3)
    endif()

    if (0)
        # Boost headers _not_ unpacked to within `include/libint2/`
        if (NOT TARGET Boost::boost)
            find_dependency(Boost 1.57)
        endif()
    endif()

    get_property(_loc TARGET ${PN}::cxx PROPERTY LOCATION)
    set(${PN}_LIBRARY ${_loc})
    get_property(_ill TARGET ${PN}::cxx PROPERTY INTERFACE_LINK_LIBRARIES)
    set(${PN}_LIBRARIES ${_ill})

    get_property(_id TARGET ${PN}::cxx PROPERTY INCLUDE_DIRECTORIES)
    set(${PN}_INCLUDE_DIR ${_id})
    get_property(_iid TARGET ${PN}::cxx PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    set(${PN}_INCLUDE_DIRS ${_iid})

    #message("Libint2::cxx")
    #message("loc ${_loc}")
    #message("ill ${_ill}")
    #message("id  ${_id}")
    #message("iid ${_iid}")
endif()
