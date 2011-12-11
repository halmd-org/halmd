#
# Define dynamic Boost Test linking flag if needed.
#
# We support both the static and dynamic library variants.
#
# http://www.boost.org/doc/libs/1_44_0/libs/test/doc/html/utf/user-guide/usage-variants.html
#
get_filename_component(
  suffix "${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}" EXT
)
if(suffix STREQUAL "${CMAKE_SHARED_LIBRARY_SUFFIX}")
  get_directory_property(
    flags COMPILE_DEFINITIONS
  )
  list(APPEND flags
    "BOOST_TEST_DYN_LINK"
  )
  set_directory_properties(PROPERTIES
    COMPILE_DEFINITIONS "${flags}"
  )
endif(suffix STREQUAL "${CMAKE_SHARED_LIBRARY_SUFFIX}")

# Link all tests against the test_tools_ctest library, which prints
# CTEST_FULL_OUTPUT to avoid ctest truncation of the test output.
set(HALMD_TEST_LIBRARIES
  test_tools_ctest
  pthread
  rt
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

set(HALMD_TEST_CMAKE_BINARY_DIR "@CMAKE_BINARY_DIR@/test/cmake")

# path to HALMD C++ and Lua modules
#
# These may be used with the ENVIRONMENT test property if needed:
#
# set_property(TEST <test>
#   PROPERTY ENVIRONMENT "LUA_PATH=${HALMD_LUA_PATH}" "LUA_CPATH=${HALMD_LUA_CPATH}"
# )
#
set(HALMD_LUA_PATH
  "${CMAKE_BINARY_DIR}/lua/?.lua\;${CMAKE_BINARY_DIR}/lua/?/init.lua"
)
set(HALMD_LUA_CPATH
  "${CMAKE_BINARY_DIR}/halmd/?${CMAKE_SHARED_LIBRARY_SUFFIX}"
)

##
# HALMD test suite
#
add_subdirectory(cmake)
add_subdirectory(integration)
add_subdirectory(lua)
add_subdirectory(performance)
add_subdirectory(tools)
add_subdirectory(unit)
