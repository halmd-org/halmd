#
# module all
#
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
  halmd_add_gpu_test(unit/mdsim/particle_groups/all/gpu
    test_unit_mdsim_particle_groups_all --run_test=gpu --log_level=test_suite
  )
endif()
add_test(unit/mdsim/particle_groups/all/host
  test_unit_mdsim_particle_groups_all --run_test=host --log_level=test_suite
)

#
# module id_range
#
add_executable(test_unit_mdsim_particle_groups_id_range
  id_range.cpp
)
target_link_libraries(test_unit_mdsim_particle_groups_id_range
  halmd_mdsim_host
  halmd_mdsim_host_particle_groups
  ${HALMD_TEST_LIBRARIES}
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_particle_groups_id_range
    halmd_mdsim_gpu_particle_groups
    halmd_mdsim_gpu
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
  halmd_add_gpu_test(unit/mdsim/particle_groups/id_range/gpu
    test_unit_mdsim_particle_groups_id_range --run_test=gpu --log_level=test_suite
  )
endif()
add_test(unit/mdsim/particle_groups/id_range/host
  test_unit_mdsim_particle_groups_id_range --run_test=host --log_level=test_suite
)

#
# module region
#
if(HALMD_WITH_GPU)
  cuda_add_executable(test_unit_mdsim_particle_groups_region
    region.cpp
    region_kernel.cu
  )
else()
  add_executable(test_unit_mdsim_particle_groups_region
    region.cpp
  )
endif()
target_link_libraries(test_unit_mdsim_particle_groups_region
  halmd_mdsim_host_particle_groups
  halmd_mdsim
  halmd_mdsim_host
  halmd_utility
  ${HALMD_TEST_LIBRARIES}
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_particle_groups_region
    halmd_mdsim_gpu_particle_groups
    halmd_mdsim_gpu
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
endif()

add_test(unit/mdsim/particle_groups/region/host/2d
  test_unit_mdsim_particle_groups_region --run_test=host/two --log_level=test_suite
)
add_test(unit/mdsim/particle_groups/region/host/3d
  test_unit_mdsim_particle_groups_region --run_test=host/three --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/mdsim/particle_groups/region/gpu/2d
    test_unit_mdsim_particle_groups_region --run_test=gpu/two --log_level=test_suite
  )
  add_test(unit/mdsim/particle_groups/region/gpu/3d
    test_unit_mdsim_particle_groups_region --run_test=gpu/three --log_level=test_suite
  )
endif()
