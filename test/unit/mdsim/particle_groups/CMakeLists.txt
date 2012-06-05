add_executable(test_unit_mdsim_particle_groups_all
  all.cpp
)
target_link_libraries(test_unit_mdsim_particle_groups_all
  halmd_mdsim_host_particle_groups
  halmd_mdsim_host
  ${HALMD_TEST_LIBRARIES}
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_particle_groups_all
    halmd_mdsim_gpu_particle_groups
    halmd_mdsim_gpu
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
  add_test(unit/mdsim/particle_groups/all/gpu
    test_unit_mdsim_particle_groups_all --run_test=gpu --log_level=test_suite
  )
endif()
add_test(unit/mdsim/particle_groups/all/host
  test_unit_mdsim_particle_groups_all --run_test=host --log_level=test_suite
)

add_test(unit/mdsim/particle_groups/from_range/host
  test_unit_mdsim_particle_groups_from_range --run_test=host --log_level=test_suite
)
add_executable(test_unit_mdsim_particle_groups_from_range
  from_range.cpp
)
target_link_libraries(test_unit_mdsim_particle_groups_from_range
  halmd_mdsim_host_particle_groups
  halmd_mdsim_host
  ${HALMD_TEST_LIBRARIES}
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_particle_groups_from_range
    halmd_mdsim_gpu_particle_groups
    halmd_mdsim_gpu
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
  add_test(unit/mdsim/particle_groups/from_range/gpu
    test_unit_mdsim_particle_groups_from_range --run_test=gpu --log_level=test_suite
  )
endif()
add_test(unit/mdsim/particle_groups/from_range/host
  test_unit_mdsim_particle_groups_from_range --run_test=host --log_level=test_suite
)
