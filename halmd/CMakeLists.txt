##
# Add HALMD library target
#
# This macro is a wrapper around add_library, which further adds the
# library target name to a global property list of HALMD libraries.
# The HALMD executable is linked against these libraries.
#
macro(halmd_add_library TARGET)
  set_property(GLOBAL APPEND PROPERTY HALMD_LIBRARIES ${TARGET})
  if(HALMD_WITH_GPU)
    cuda_add_library(${TARGET} ${ARGN})
  else()
    add_library(${TARGET} ${ARGN})
  endif()
endmacro()

macro(halmd_add_executable)
  if(HALMD_WITH_GPU)
    cuda_add_executable(${ARGN})
  else()
    add_executable(${ARGN})
  endif()
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
# Add optional HALMD potential
#
# This macro automatically adds the specified potential module if enabled via
# the CMake cache variable HALMD_POTENTIALS, see cmake/potentials.cmake.
#
# The first argument is the library name, the second the type of the potential
# (e.g., pair), and the third one the name of the potential, the remaining
# arguments are passed to add_library()
#
macro(halmd_add_potential LIBRARY TYPE POTENTIAL)
  if(HALMD_WITH_${TYPE}_${POTENTIAL})
    halmd_add_library(${LIBRARY} ${ARGN})
    halmd_add_modules("lib${LIBRARY}")
  endif()
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

halmd_add_library(halmd_script
  script.cpp
)

##
# Generate Lua binding function calls from module target sources
#
macro(halmd_bundle_modules INPUT OUTPUT LIBRARY)
  set(LIBRARY "${LIBRARY}")
  set(MODULE)
  file(STRINGS "${INPUT}" lines REGEX "@MODULE@" LIMIT_COUNT 1)
  foreach(line ${lines})
    if(line MATCHES "^(.*)@MODULE@(.*)$")
      string(REPLACE ";" "${CMAKE_MATCH_2}\n${CMAKE_MATCH_1}" MODULE "${ARGN}")
    endif()
  endforeach()
  configure_file("${INPUT}" "${OUTPUT}" @ONLY)
endmacro()

# Store list of all HALMD C++ modules into variable
get_property(HALMD_MODULES GLOBAL PROPERTY HALMD_MODULES)
# Store list of all HALMD C++ libraries into variable
get_property(HALMD_LIBRARIES GLOBAL PROPERTY HALMD_LIBRARIES)

# Create header file with luaopen_halmd_<module> declarations
halmd_bundle_modules(modules.hpp.in modules.hpp halmd ${HALMD_MODULES})
# Create source file with luaopen_halmd_<module> calls
halmd_bundle_modules(modules.cpp.in modules.cpp halmd ${HALMD_MODULES})

# Generate configuration header
configure_file(config.hpp.in config.hpp @ONLY)
# Substitute version info in header file
configure_file(version.h.in version.h ESCAPE_QUOTES @ONLY)

##
# HALMD executable
#
halmd_add_executable(halmd
  main.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/modules.cpp"
)
target_link_libraries(halmd
  # with newer version of GNU ld, indirect linking through intermediate
  # libraries does not work anymore, so we have to explicitly link halmd
  # against common libraries
  # http://fedoraproject.org/wiki/UnderstandingDSOLinkChange
  ${HALMD_LIBRARIES}
  ${HALMD_COMMON_LIBRARIES}
)
install(TARGETS halmd
  RUNTIME DESTINATION bin
)
