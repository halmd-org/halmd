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

# Test: colloidal particle in a harmonic trap
add_executable(test_integration_validation_trapped_colloid
  trapped_colloid.cpp
)
target_link_libraries(test_integration_validation_trapped_colloid
  ${HALMD_TEST_LIBRARIES}
)

configure_file(
  trapped_colloid.cmake.in
  test_integration_validation_trapped_colloid.cmake
  @ONLY
)
halmd_add_test(integration/validation/trapped_colloid/host/3d
  ${CMAKE_COMMAND} -DDIMENSION=3 -DUSE_GPU=0 -P test_integration_validation_trapped_colloid.cmake
)
halmd_add_test(integration/validation/trapped_colloid/host/2d
  ${CMAKE_COMMAND} -DDIMENSION=2 -DUSE_GPU=0 -P test_integration_validation_trapped_colloid.cmake
)
if(HALMD_WITH_GPU)
  halmd_add_gpu_test(NO_MEMCHECK integration/validation/trapped_colloid/gpu/3d
    ${CMAKE_COMMAND} -DDIMENSION=3 -DUSE_GPU=1 -P test_integration_validation_trapped_colloid.cmake
  )
  halmd_add_gpu_test(NO_MEMCHECK integration/validation/trapped_colloid/gpu/2d
    ${CMAKE_COMMAND} -DDIMENSION=2 -DUSE_GPU=1 -P test_integration_validation_trapped_colloid.cmake
  )
endif()
