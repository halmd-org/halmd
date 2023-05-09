add_executable(test_unit_mdsim_orientations_uniform
  uniform.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_orientations_uniform
    halmd_mdsim_gpu_orientations
    halmd_mdsim_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_orientations_uniform
  halmd_mdsim_host_orientations
  halmd_mdsim_host
  halmd_mdsim
  halmd_observables
  halmd_random_host
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/orientations/uniform/host/2d
  test_unit_mdsim_orientations_uniform --run_test=uniform_host_2d --log_level=test_suite
)
add_test(unit/mdsim/orientations/uniform/host/3d
  test_unit_mdsim_orientations_uniform --run_test=uniform_host_3d --log_level=test_suite
)
if(HALMD_WITH_GPU)
  if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
    halmd_add_gpu_test(unit/mdsim/orientations/uniform/gpu/float/2d
      test_unit_mdsim_orientations_uniform --run_test=uniform_gpu_float_2d --log_level=test_suite
    )
    halmd_add_gpu_test(unit/mdsim/orientations/uniform/gpu/float/3d
      test_unit_mdsim_orientations_uniform --run_test=uniform_gpu_float_3d --log_level=test_suite
    )
  endif()
  if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
    halmd_add_gpu_test(unit/mdsim/orientations/uniform/gpu/dsfloat/2d
      test_unit_mdsim_orientations_uniform --run_test=uniform_gpu_dsfloat_2d --log_level=test_suite
    )
    halmd_add_gpu_test(unit/mdsim/orientations/uniform/gpu/dsfloat/3d
      test_unit_mdsim_orientations_uniform --run_test=uniform_gpu_dsfloat_3d --log_level=test_suite
    )
  endif()
endif()
