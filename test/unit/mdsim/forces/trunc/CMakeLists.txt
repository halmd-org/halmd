if(HALMD_WITH_pair_lennard_jones)
  add_executable(test_unit_mdsim_forces_trunc_smooth_r4
          smooth_r4.cpp
  )
  if(HALMD_WITH_GPU)
    target_link_libraries(test_unit_mdsim_forces_trunc_smooth_r4
      halmd_mdsim_gpu_potentials_pair_lennard_jones
      halmd_mdsim_gpu
      halmd_mdsim
      halmd_algorithm_gpu
      halmd_utility_gpu
    )
  endif()
  target_link_libraries(test_unit_mdsim_forces_trunc_smooth_r4
    halmd_mdsim_host_potentials_pair_lennard_jones
    ${HALMD_TEST_LIBRARIES}
  )
  add_test(unit/mdsim/forces/trunc/smooth_r4/host
    test_unit_mdsim_forces_trunc_smooth_r4 --run_test=host/smooth_r4 --log_level=test_suite
  )
  if(HALMD_WITH_GPU)
    halmd_add_gpu_test(unit/mdsim/forces/trunc/smooth_r4/gpu
      test_unit_mdsim_forces_trunc_smooth_r4 --run_test=gpu/smooth_r4 --log_level=test_suite
    )
  endif()
endif()
