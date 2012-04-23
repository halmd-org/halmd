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

# Split HALMD C++ modules into GPU and non-GPU modules
set(HALMD_GPU_MODULES)
foreach(MODULE ${HALMD_MODULES})
  if(MODULE MATCHES "gpu")
    list(APPEND HALMD_GPU_MODULES ${MODULE})
    list(REMOVE_ITEM HALMD_MODULES ${MODULE})
  endif()
endforeach()

# Split HALMD C++ libraries into GPU and non-GPU libraries
set(HALMD_GPU_LIBRARIES)
foreach(LIBRARY ${HALMD_LIBRARIES})
  if(LIBRARY MATCHES "gpu")
    list(APPEND HALMD_GPU_LIBRARIES ${LIBRARY})
    list(REMOVE_ITEM HALMD_LIBRARIES ${LIBRARY})
  endif()
endforeach()

# Generate configuration header
configure_file(config.hpp.in config.hpp @ONLY)

# Substitute version info in header file
configure_file(version.h.in version.h ESCAPE_QUOTES @ONLY)

halmd_bundle_modules(modules.hpp.in modules.hpp halmd ${HALMD_MODULES} ${HALMD_GPU_MODULES})

##
# HALMD executable with non-GPU modules
#

halmd_bundle_modules(modules.cpp.in modules_halmd.cpp halmd ${HALMD_MODULES})

add_executable(halmd
  main.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/modules_halmd.cpp"
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

##
# HALMD plugin with GPU modules
#

if(WITH_CUDA)
  halmd_bundle_modules(modules.cpp.in modules_halmd_gpu.cpp halmd_gpu ${HALMD_GPU_MODULES})

  add_library(halmd_gpu MODULE
    "${CMAKE_CURRENT_BINARY_DIR}/modules_halmd_gpu.cpp"
  )
  set_target_properties(halmd_gpu PROPERTIES
    PREFIX ""
    OUTPUT_NAME "gpu"
    # http://www.cmake.org/Wiki/CMake_FAQ#Why_are_libraries_linked_to_my_shared_library_included_when_something_links_to_it.3F
    LINK_INTERFACE_LIBRARIES ""
  )
  if(CMAKE_CXX_PLATFORM_ID STREQUAL "Linux")
    set_target_properties(halmd_gpu PROPERTIES LINK_FLAGS "-Wl,--as-needed")
  endif()
  target_link_libraries(halmd_gpu
    ${HALMD_GPU_LIBRARIES}
    ${HALMD_COMMON_LIBRARIES}
  )
  install(TARGETS halmd_gpu
    LIBRARY DESTINATION lib/halmd
  )
endif()
