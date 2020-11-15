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
  halmd_io
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  ${HALMD_COMMON_LIBRARIES}
)
# make sure to link against halmd_utility_gpu, required for halmd::device d'tor
if(HALMD_WITH_GPU)
  list(APPEND HALMD_TEST_LIBRARIES
    halmd_utility_gpu
  )
endif(HALMD_WITH_GPU)

set(HALMD_TEST_CMAKE_BINARY_DIR "${CMAKE_BINARY_DIR}/test/cmake")

# path to HALMD executable
set(HALMD_EXECUTABLE "${CMAKE_BINARY_DIR}/halmd/halmd")

# path to HALMD C++ and Lua modules
#
# These may be used with the ENVIRONMENT test property if needed:
#
# set_property(TEST <test>
#   PROPERTY ENVIRONMENT "LUA_PATH=${HALMD_LUA_PATH}" "LUA_CPATH=${HALMD_LUA_CPATH}"
# )
#
set(HALMD_LUA_PATH
  "${CMAKE_BINARY_DIR}/lua/?.lua\;${CMAKE_BINARY_DIR}/lua/?/init.lua\;${CMAKE_CURRENT_SOURCE_DIR}/tools/?.lua"
)
set(HALMD_LUA_CPATH
  "${CMAKE_BINARY_DIR}/halmd/?${CMAKE_SHARED_LIBRARY_SUFFIX}"
)

if(HALMD_WITH_GPU)
  option(HALMD_TEST_CUDA_MEMCHECK "Set to OFF to disable CUDA memory checks while running the tests." ON)

  if(HALMD_TEST_CUDA_MEMCHECK)
    find_program(CUDA_MEMCHECK cuda-memcheck)
    if(CUDA_MEMCHECK)
      set(CUDA_MEMCHECK_COMMAND ${CUDA_MEMCHECK} --error-exitcode 1)
    else()
      message(WARNING "cuda-memcheck was not found. No memory checks will be performed.")
      set(CUDA_MEMCHECK_COMMAND)
    endif()
  else()
    set(CUDA_MEMCHECK_COMMAND)
  endif()

  macro(halmd_add_gpu_test NAME)
    if ("${NAME}" STREQUAL "NO_MEMCHECK")
      add_test(${ARGN})
    else()
      add_test(${NAME} ${CUDA_MEMCHECK_COMMAND} ${ARGN})
    endif()
  endmacro()
endif()

##
# HALMD test suite
#
add_subdirectory(cmake)
add_subdirectory(integration)
add_subdirectory(lua)
add_subdirectory(performance)
add_subdirectory(tools)
add_subdirectory(unit)
