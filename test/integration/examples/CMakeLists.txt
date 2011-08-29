configure_file(
  thermostat.cmake.in
  test_integration_examples_thermostat.cmake
  @ONLY
)
add_test(integration/examples/thermostat/host/2d
  @CMAKE_COMMAND@ -DDIMENSION=2 -DUSE_HOST=1 -P test_integration_examples_thermostat.cmake
)
add_test(integration/examples/thermostat/host/3d
  @CMAKE_COMMAND@ -DDIMENSION=3 -DUSE_HOST=1 -P test_integration_examples_thermostat.cmake
)
if(WITH_CUDA)
  add_test(integration/examples/thermostat/gpu/2d
    @CMAKE_COMMAND@ -DDIMENSION=2 -P test_integration_examples_thermostat.cmake
  )
  add_test(integration/examples/thermostat/gpu/3d
    @CMAKE_COMMAND@ -DDIMENSION=3 -P test_integration_examples_thermostat.cmake
  )
endif()
