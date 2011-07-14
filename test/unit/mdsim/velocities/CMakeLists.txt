add_executable(test_unit_mdsim_velocities_boltzmann
  boltzmann.cpp
)
target_link_libraries(test_unit_mdsim_velocities_boltzmann
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/velocities/boltzmann/host/2d
  test_unit_mdsim_velocities_boltzmann --run_test=boltzmann_host_2d --log_level=test_suite
)
add_test(unit/mdsim/velocities/boltzmann/host/3d
  test_unit_mdsim_velocities_boltzmann --run_test=boltzmann_host_3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/velocities/boltzmann/gpu/2d
    test_unit_mdsim_velocities_boltzmann --run_test=boltzmann_gpu_2d --log_level=test_suite
  )
  add_test(unit/mdsim/velocities/boltzmann/gpu/3d
    test_unit_mdsim_velocities_boltzmann --run_test=boltzmann_gpu_3d --log_level=test_suite
  )
endif(WITH_CUDA)
