#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "vb" for configuration "Release"
set_property(TARGET vb APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(vb PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_RELEASE "CXX"
  IMPORTED_LOCATION_RELEASE "/usr/local/vb/lib/libvb.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS vb )
list(APPEND _IMPORT_CHECK_FILES_FOR_vb "/usr/local/vb/lib/libvb.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
