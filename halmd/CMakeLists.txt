##
# Add HALMD library target
#
# This macro is a wrapper around add_library, which further adds the
# library target name to a global property list of HALMD libraries.
# The HALMD executable is linked against these libraries.
#
macro(halmd_add_library TARGET)
  set_property(GLOBAL APPEND PROPERTY HALMD_LIBRARIES ${TARGET})
  add_library(${TARGET} ${ARGN})
endmacro()

##
# Add HALMD modules
#
# This macro adds all arguments to a global property list of HALMD modules.
# This list is used to generate calls to luaopen_<module> functions, which
# exports a C++ module into the Lua interpreter.
#
# Base classes must be ordered *before* derived classes, otherwise
# Luabind will throw an assertion error or cause a segmentation fault.
#
macro(halmd_add_modules)
  set_property(GLOBAL APPEND PROPERTY HALMD_MODULES ${ARGN})
endmacro()

##
# HALMD C++ modules
#
add_subdirectory(algorithm)
add_subdirectory(io)
add_subdirectory(mdsim)
add_subdirectory(observables)
add_subdirectory(random)
add_subdirectory(utility)

##
# Generate Lua binding function calls from module target sources
#
get_property(HALMD_MODULES GLOBAL PROPERTY HALMD_MODULES)

foreach(source modules.hpp modules.cpp)
  file(STRINGS "${source}.in" lines REGEX "@MODULE@" LIMIT_COUNT 1)
  foreach(line ${lines})
    if(line MATCHES "^(.*)@MODULE@(.*)$")
      string(REPLACE ";" "${CMAKE_MATCH_2}\n${CMAKE_MATCH_1}" MODULE "${HALMD_MODULES}")
    endif()
  endforeach()
  set(LIBRARY libhalmd)
  configure_file("${source}.in" "${source}" @ONLY)
  set(MODULE)
endforeach()

configure_file(version.h.in version.h ESCAPE_QUOTES @ONLY)

##
# Add HALMD executable target
#
get_property(HALMD_LIBRARIES GLOBAL PROPERTY HALMD_LIBRARIES)

add_library(libhalmd SHARED
  "${CMAKE_CURRENT_BINARY_DIR}/modules.cpp"
  script.cpp
)
set_target_properties(libhalmd PROPERTIES
  PREFIX ""
  # http://www.cmake.org/Wiki/CMake_FAQ#Why_are_libraries_linked_to_my_shared_library_included_when_something_links_to_it.3F
  LINK_INTERFACE_LIBRARIES ""
)
target_link_libraries(libhalmd
  ${HALMD_LIBRARIES}
  ${HALMD_COMMON_LIBRARIES}
)
install(TARGETS libhalmd
  LIBRARY DESTINATION lib
)

add_executable(halmd
  main.cpp
)
target_link_libraries(halmd
  libhalmd
  # needed for halmd::host_name
  halmd_utility
  # with newer version of GNU ld, indirect linking through intermediate
  # libraries does not work anymore, so we have to explicitly link halmd
  # against common libraries
  # http://fedoraproject.org/wiki/UnderstandingDSOLinkChange
  ${HALMD_COMMON_LIBRARIES}
)
install(TARGETS halmd
  RUNTIME DESTINATION bin
)

# On Linux, add --Wl,--as-needed to the link flags of the halmd executable.
# This makes distribution packagers happy, as the halmd binary will only
# reference the actually used libraries, which allows for minimal package
# dependencies of the halmd package. The HALMD_COMMON_LIBRARIES variable
# contains more libraries than needed in most cases, to support linking
# Boost and HDF5 libraries statically.
#
# Note that -Wl,--as-needed must not be set globally, the unit tests
# have to be build _without_ this flag. Otherwise fixtures in linked
# shared test libraries will be silently discarded by the linker.
#
if(CMAKE_CXX_PLATFORM_ID STREQUAL "Linux")
  set_target_properties(halmd PROPERTIES LINK_FLAGS "-Wl,--as-needed")
endif()
