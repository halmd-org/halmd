add_executable(test_unit_mdsim_smoothers_localr4
  localr4.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_smoothers_localr4
    halmd_mdsim_gpu_forces
    halmd_mdsim_gpu_potentials
    halmd_mdsim_gpu
    halmd_mdsim
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_smoothers_localr4
  halmd_mdsim_host_potentials
  halmd_mdsim_smoothers
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/smoothers/localr4/host
  test_unit_mdsim_smoothers_localr4 --run_test=host/localr4 --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/mdsim/smoothers/localr4/gpu
    test_unit_mdsim_smoothers_localr4 --run_test=gpu/localr4 --log_level=test_suite
  )
endif()
