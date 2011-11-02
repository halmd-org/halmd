add_executable(test_unit_mdsim_integrators_verlet
  verlet.cpp
)
target_link_libraries(test_unit_mdsim_integrators_verlet
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/integrators/verlet/host/2d
  test_unit_mdsim_integrators_verlet --run_test=ideal_gas_host_2d --log_level=test_suite
)
add_test(unit/mdsim/integrators/verlet/host/3d
  test_unit_mdsim_integrators_verlet --run_test=ideal_gas_host_3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/integrators/verlet/gpu/2d
    test_unit_mdsim_integrators_verlet --run_test=ideal_gas_gpu_2d --log_level=test_suite
  )
  add_test(unit/mdsim/integrators/verlet/gpu/3d
    test_unit_mdsim_integrators_verlet --run_test=ideal_gas_gpu_3d --log_level=test_suite
  )
endif(WITH_CUDA)

add_executable(test_unit_mdsim_integrators_euler
  euler.cpp
)
target_link_libraries(test_unit_mdsim_integrators_euler
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/integrators/euler/host/2d/linear
  test_unit_mdsim_integrators_euler --run_test=euler_host_2d_linear --log_level=test_suite
)
add_test(unit/mdsim/integrators/euler/host/2d/overdamped
  test_unit_mdsim_integrators_euler --run_test=euler_host_2d_overdamped --log_level=test_suite
)
add_test(unit/mdsim/integrators/euler/host/3d/linear
  test_unit_mdsim_integrators_euler --run_test=euler_host_3d_linear --log_level=test_suite
)
add_test(unit/mdsim/integrators/euler/host/3d/overdamped
  test_unit_mdsim_integrators_euler --run_test=euler_host_3d_overdamped --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/integrators/euler/gpu/2d/linear
    test_unit_mdsim_integrators_euler --run_test=euler_gpu_2d_linear --log_level=test_suite
  )
  add_test(unit/mdsim/integrators/euler/gpu/2d/overdamped
    test_unit_mdsim_integrators_euler --run_test=euler_gpu_2d_overdamped --log_level=test_suite
  )
  add_test(unit/mdsim/integrators/euler/gpu/3d/linear
    test_unit_mdsim_integrators_euler --run_test=euler_gpu_3d_linear --log_level=test_suite
  )
  add_test(unit/mdsim/integrators/euler/gpu/3d/overdamped
    test_unit_mdsim_integrators_euler --run_test=euler_gpu_3d_overdamped --log_level=test_suite
  )
endif(WITH_CUDA)

add_executable(test_unit_mdsim_integrators_verlet_nvt_andersen
  verlet_nvt_andersen.cpp
)
target_link_libraries(test_unit_mdsim_integrators_verlet_nvt_andersen
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/integrators/verlet_nvt_andersen/host/2d
  test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen_host_2d --log_level=test_suite
)
add_test(unit/mdsim/integrators/verlet_nvt_andersen/host/3d
  test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen_host_3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/integrators/verlet_nvt_andersen/gpu/2d
    test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen_gpu_2d --log_level=test_suite
  )
  add_test(unit/mdsim/integrators/verlet_nvt_andersen/gpu/3d
    test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen_gpu_3d --log_level=test_suite
  )
endif(WITH_CUDA)
