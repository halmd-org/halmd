add_executable(test_unit_mdsim_forces_interpolation_cubic_hermite
  cubic_hermite.cpp
)

add_executable(test_unit_mdsim_forces_interpolation_linear
  linear.cpp
)

if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_forces_interpolation_cubic_hermite
    halmd_mdsim_gpu
    halmd_mdsim
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_forces_interpolation_cubic_hermite
  halmd_mdsim_forces_interpolation
  halmd_mdsim
  ${HALMD_TEST_LIBRARIES}
)

target_link_libraries(test_unit_mdsim_forces_interpolation_linear
  halmd_mdsim_forces_interpolation
  halmd_mdsim
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/forces/interpolation/cubic_hermite/host
  test_unit_mdsim_forces_interpolation_cubic_hermite --run_test=host/cubic_hermite --log_level=test_suite
)

add_test(unit/mdsim/forces/interpolation/linear/host
  test_unit_mdsim_forces_interpolation_linear --run_test=host/linear --log_level=test_suite
)
if(HALMD_WITH_GPU)
#  add_test(unit/mdsim/forces/interpolation/cubic_hermite/gpu
#    test_unit_mdsim_forces_interpolation_cubic_hermite --run_test=gpu/cubic_hermite --log_level=test_suite
#  )
endif()
