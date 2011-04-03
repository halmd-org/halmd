##
# Add HALMD module targets
#

include(halmd_module)

halmd_add_module(halmd_runner
  runner.cpp
)

add_subdirectory(algorithm)
add_subdirectory(io)
add_subdirectory(mdsim)
add_subdirectory(observables)
add_subdirectory(random)
add_subdirectory(utility)

##
# Generate Lua binding function calls from list of modules
#

get_property(HALMD_MODULES GLOBAL PROPERTY HALMD_MODULES)

foreach(source modules.hpp modules.cpp)
  file(STRINGS "${source}.in" lines REGEX "@MODULE_NAME@" LIMIT_COUNT 1)
  foreach(line ${lines})
    if(line MATCHES "^(.*)@MODULE_NAME@(.*)$")
      string(REPLACE ";" "${CMAKE_MATCH_2}\n${CMAKE_MATCH_1}" MODULE_NAME "${HALMD_MODULES}")
    endif()
  endforeach()
  configure_file("${source}.in" "${source}" @ONLY)
  set(MODULE_NAME)
endforeach()

##
# Add HALMD executable target
#

configure_file(version.h.in version.h ESCAPE_QUOTES @ONLY)

add_executable(halmd
  main.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/modules.cpp"
  script.cpp
)
target_link_libraries(halmd
  ${HALMD_MODULES}
  halmd_utility_hostname
  ${HALMD_COMMON_LIBRARIES}
)
install(TARGETS halmd
  RUNTIME DESTINATION bin
)
