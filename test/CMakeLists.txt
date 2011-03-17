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
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)

#
# Include unit tests
#

add_subdirectory(algorithm)
add_subdirectory(h5xx)
add_subdirectory(io)
add_subdirectory(mdsim)
add_subdirectory(numeric)
add_subdirectory(observables)
add_subdirectory(performance)
add_subdirectory(random)
add_subdirectory(tools)
add_subdirectory(utility)

if(WITH_CUDA)
  add_subdirectory(cuda_wrapper)
endif(WITH_CUDA)
