#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Libint2::cxx" for configuration "Release"
set_property(TARGET Libint2::cxx APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Libint2::cxx PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libint2.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS Libint2::cxx )
list(APPEND _IMPORT_CHECK_FILES_FOR_Libint2::cxx "${_IMPORT_PREFIX}/lib/libint2.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
