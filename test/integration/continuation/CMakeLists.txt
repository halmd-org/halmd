##
# HALMD trajectory continuation tests
#
configure_file(
  lennard_jones.cmake.in
  test_integration_continuation_lennard_jones.cmake
  @ONLY
)
add_test(integration/continuation/lennard_jones
  test_integration_continuation_lennard_jones.cmake
)
