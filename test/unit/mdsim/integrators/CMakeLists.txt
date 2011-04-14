# module verlet_nvt_andersen
add_executable(test_unit_mdsim_integrators_verlet_nvt_andersen
  verlet_nvt_andersen.cpp
)
target_link_libraries(test_unit_mdsim_integrators_verlet_nvt_andersen
  halmd_mdsim_host_forces
  halmd_mdsim_host_integrators
  halmd_mdsim_host_positions
  halmd_mdsim_host_velocities
  halmd_mdsim_host
  halmd_observables_host
  halmd_random_host
)
if(WITH_CUDA)
  target_link_libraries(test_unit_mdsim_integrators_verlet_nvt_andersen
    halmd_mdsim_gpu_forces
    halmd_mdsim_gpu_integrators
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu_velocities
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif(WITH_CUDA)
target_link_libraries(test_unit_mdsim_integrators_verlet_nvt_andersen
  halmd_io
  halmd_mdsim
  halmd_random
  halmd_observables
  halmd_utility
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/integrators/verlet_nvt_andersen/host/2d
  test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen/host/2d --log_level=test_suite
)
add_test(unit/mdsim/integrators/verlet_nvt_andersen/host/3d
  test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen/host/3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/integrators/verlet_nvt_andersen/gpu/2d
    test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen/gpu/2d --log_level=test_suite
  )
  add_test(unit/mdsim/integrators/verlet_nvt_andersen/gpu/3d
    test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen/gpu/3d --log_level=test_suite
  )
endif(WITH_CUDA)


# module verlet_nvt_hoover
add_executable(test_unit_mdsim_integrators_verlet_nvt_hoover
  verlet_nvt_hoover.cpp
)
target_link_libraries(test_unit_mdsim_integrators_verlet_nvt_hoover
  halmd_mdsim_host_forces
  halmd_mdsim_host_integrators
  halmd_mdsim_host_positions
  halmd_mdsim_host_velocities
  halmd_mdsim_host
  halmd_observables_host
  halmd_random_host
)
if(WITH_CUDA)
  target_link_libraries(test_unit_mdsim_integrators_verlet_nvt_hoover
    halmd_mdsim_gpu_forces
    halmd_mdsim_gpu_integrators
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu_velocities
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif(WITH_CUDA)
target_link_libraries(test_unit_mdsim_integrators_verlet_nvt_hoover
  halmd_io
  halmd_mdsim
  halmd_random
  halmd_observables
  halmd_utility
  ${HALMD_COMMON_LIBRARIES}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(unit/mdsim/integrators/verlet_nvt_hoover/host/2d
  test_unit_mdsim_integrators_verlet_nvt_hoover --run_test=verlet_nvt_hoover/host/2d --log_level=test_suite
)
add_test(unit/mdsim/integrators/verlet_nvt_hoover/host/3d
  test_unit_mdsim_integrators_verlet_nvt_hoover --run_test=verlet_nvt_hoover/host/3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/integrators/verlet_nvt_hoover/gpu/2d
    test_unit_mdsim_integrators_verlet_nvt_hoover --run_test=verlet_nvt_hoover/gpu/2d --log_level=test_suite
  )
  add_test(unit/mdsim/integrators/verlet_nvt_hoover/gpu/3d
    test_unit_mdsim_integrators_verlet_nvt_hoover --run_test=verlet_nvt_hoover/gpu/3d --log_level=test_suite
  )
endif(WITH_CUDA)
