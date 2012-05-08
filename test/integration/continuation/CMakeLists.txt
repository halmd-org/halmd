##
# HALMD trajectory continuation tests
#
configure_file(
  lennard_jones.cmake.in
  test_integration_continuation_lennard_jones.cmake
  @ONLY
)
add_test(integration/continuation/lennard_jones/host/2d
  @CMAKE_COMMAND@ -DDIMENSION=2 -DUSE_HOST=1 -P test_integration_continuation_lennard_jones.cmake
)
add_test(integration/continuation/lennard_jones/host/3d
  @CMAKE_COMMAND@ -DDIMENSION=3 -DUSE_HOST=1 -P test_integration_continuation_lennard_jones.cmake
)
if(HALMD_WITH_GPU)
  add_test(integration/continuation/lennard_jones/gpu/2d
    @CMAKE_COMMAND@ -DDIMENSION=2 -P test_integration_continuation_lennard_jones.cmake
  )
  add_test(integration/continuation/lennard_jones/gpu/3d
    @CMAKE_COMMAND@ -DDIMENSION=3 -P test_integration_continuation_lennard_jones.cmake
  )
endif()
configure_file(
  lennard_jones_simple.cmake.in
  test_integration_continuation_lennard_jones_simple.cmake
  @ONLY
)
add_test(integration/continuation/lennard_jones_simple/host/2d
  @CMAKE_COMMAND@ -DDIMENSION=2 -DUSE_HOST=1 -P test_integration_continuation_lennard_jones_simple.cmake
)
add_test(integration/continuation/lennard_jones_simple/host/3d
  @CMAKE_COMMAND@ -DDIMENSION=3 -DUSE_HOST=1 -P test_integration_continuation_lennard_jones_simple.cmake
)
if(HALMD_WITH_GPU)
  add_test(integration/continuation/lennard_jones_simple/gpu/2d
    @CMAKE_COMMAND@ -DDIMENSION=2 -P test_integration_continuation_lennard_jones_simple.cmake
  )
  add_test(integration/continuation/lennard_jones_simple/gpu/3d
    @CMAKE_COMMAND@ -DDIMENSION=3 -P test_integration_continuation_lennard_jones_simple.cmake
  )
endif()
