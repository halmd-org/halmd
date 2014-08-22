add_subdirectory(forces)
add_subdirectory(integrators)
add_subdirectory(particle_groups)
add_subdirectory(positions)
add_subdirectory(potentials)
add_subdirectory(velocities)

add_executable(test_unit_mdsim_binning
  binning.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_binning
    halmd_mdsim_gpu
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_binning
  halmd_mdsim_host
  halmd_mdsim
  halmd_utility
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/binning/host/2d
  test_unit_mdsim_binning --run_test=host/two --log_level=test_suite
)
add_test(unit/mdsim/binning/host/3d
  test_unit_mdsim_binning --run_test=host/three --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/mdsim/binning/gpu/2d
    test_unit_mdsim_binning --run_test=gpu/two --log_level=test_suite
  )
  add_test(unit/mdsim/binning/gpu/3d
    test_unit_mdsim_binning --run_test=gpu/three --log_level=test_suite
  )
endif()

# module box
if(HALMD_WITH_GPU)
  add_executable(test_unit_mdsim_box
    box.cpp
    box_kernel.cu
  )
else (HALMD_WITH_GPU)
  add_executable(test_unit_mdsim_box
    box.cpp
  )
endif()
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_box
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_box
  halmd_mdsim
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/box/construction/2d
  test_unit_mdsim_box --run_test=box_construction_2d --log_level=test_suite
)
add_test(unit/mdsim/box/construction/3d
  test_unit_mdsim_box --run_test=box_construction_3d --log_level=test_suite
)
add_test(unit/mdsim/box/periodic/host/2d
  test_unit_mdsim_box --run_test=box_periodic_host_2d --log_level=test_suite
)
add_test(unit/mdsim/box/periodic/host/3d
  test_unit_mdsim_box --run_test=box_periodic_host_3d --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/mdsim/box/periodic/gpu/2d
    test_unit_mdsim_box --run_test=box_periodic_gpu_2d --log_level=test_suite
  )
  add_test(unit/mdsim/box/periodic/gpu/3d
    test_unit_mdsim_box --run_test=box_periodic_gpu_3d --log_level=test_suite
  )
endif()

# module clock
add_executable(test_unit_mdsim_clock
  clock.cpp
)
target_link_libraries(test_unit_mdsim_clock
  halmd_mdsim
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/clock
  test_unit_mdsim_clock --log_level=test_suite
)

# module particle
add_executable(test_unit_mdsim_particle
  particle.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_particle
    halmd_mdsim_gpu
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_particle
  halmd_mdsim_host
  ${HALMD_TEST_LIBRARIES}
)

add_test(unit/mdsim/particle/host/2d
  test_unit_mdsim_particle --run_test=host/two --log_level=test_suite
)
add_test(unit/mdsim/particle/host/3d
  test_unit_mdsim_particle --run_test=host/three --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/mdsim/particle/gpu/2d
    test_unit_mdsim_particle --run_test=gpu/two --log_level=test_suite
  )
  add_test(unit/mdsim/particle/gpu/3d
    test_unit_mdsim_particle --run_test=gpu/three --log_level=test_suite
  )
endif()

# module region
if(HALMD_WITH_GPU)
  add_executable(test_unit_mdsim_region
    region.cpp
    region_kernel.cu
  )
else (HALMD_WITH_GPU)
  add_executable(test_unit_mdsim_region
    region.cpp
  )
endif()
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_region
    halmd_mdsim_gpu
    halmd_algorithm_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_region
  halmd_mdsim
  halmd_mdsim_host
  halmd_utility
  ${HALMD_TEST_LIBRARIES}
)

add_test(unit/mdsim/region/host/2d
  test_unit_mdsim_region --run_test=host/two --log_level=test_suite
)
add_test(unit/mdsim/region/host/3d
  test_unit_mdsim_region --run_test=host/three --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/mdsim/region/gpu/2d
    test_unit_mdsim_region --run_test=gpu/two --log_level=test_suite
  )
  add_test(unit/mdsim/region/gpu/3d
    test_unit_mdsim_region --run_test=gpu/three --log_level=test_suite
  )
endif()
