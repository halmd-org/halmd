##
# HALMD validation tests
#
# Thermodynamics Test
add_executable(test_integration_validation_thermodynamics
  thermodynamics.cpp
)
target_link_libraries(test_integration_validation_thermodynamics
  ${HALMD_TEST_LIBRARIES}
)

configure_file(
  thermodynamics.cmake.in
  test_integration_validation_thermodynamics.cmake
  @ONLY
)
halmd_add_test(integration/validation/thermodynamics/host/2d
  ${CMAKE_COMMAND} -DDIMENSION=2 -DUSE_GPU=0 -P test_integration_validation_thermodynamics.cmake
)
halmd_add_test(integration/validation/thermodynamics/host/3d
  ${CMAKE_COMMAND} -DDIMENSION=3 -DUSE_GPU=0 -P test_integration_validation_thermodynamics.cmake
)
if(HALMD_WITH_GPU)
  halmd_add_gpu_test(NO_MEMCHECK integration/validation/thermodynamics/gpu/2d
    ${CMAKE_COMMAND} -DDIMENSION=2 -DUSE_GPU=1 -P test_integration_validation_thermodynamics.cmake
  )
  halmd_add_gpu_test(NO_MEMCHECK integration/validation/thermodynamics/gpu/3d
    ${CMAKE_COMMAND} -DDIMENSION=3 -DUSE_GPU=1 -P test_integration_validation_thermodynamics.cmake
  )
endif()

# Ornstein-Uhlenbeck Validation Test
add_executable(test_integration_validation_ornstein_uhlenbeck
  ornstein_uhlenbeck.cpp
)
target_link_libraries(test_integration_validation_ornstein_uhlenbeck
  ${HALMD_TEST_LIBRARIES}
)

configure_file(
  ornstein_uhlenbeck.cmake.in
  test_integration_validation_ornstein_uhlenbeck.cmake
  @ONLY
)
halmd_add_test(integration/validation/ornstein_uhlenbeck/host/3d
  ${CMAKE_COMMAND} -DDIMENSION=3 -DUSE_GPU=0 -P test_integration_validation_ornstein_uhlenbeck.cmake
)
halmd_add_test(integration/validation/ornstein_uhlenbeck/host/2d
  ${CMAKE_COMMAND} -DDIMENSION=2 -DUSE_GPU=0 -P test_integration_validation_ornstein_uhlenbeck.cmake
)
if(HALMD_WITH_GPU)
  halmd_add_gpu_test(NO_MEMCHECK integration/validation/ornstein_uhlenbeck/gpu/3d
    ${CMAKE_COMMAND} -DDIMENSION=3 -DUSE_GPU=1 -P test_integration_validation_ornstein_uhlenbeck.cmake
  )
  halmd_add_gpu_test(NO_MEMCHECK integration/validation/ornstein_uhlenbeck/gpu/2d
    ${CMAKE_COMMAND} -DDIMENSION=2 -DUSE_GPU=1 -P test_integration_validation_ornstein_uhlenbeck.cmake
  )
endif()
