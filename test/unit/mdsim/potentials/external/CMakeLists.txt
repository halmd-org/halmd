if(${HALMD_WITH_external_planar_wall})
  add_executable(test_unit_mdsim_potentials_external_planar_wall
    planar_wall.cpp
  )
  if(HALMD_WITH_GPU)
    target_link_libraries(test_unit_mdsim_potentials_external_planar_wall
      halmd_mdsim_gpu_potentials_external_planar_wall
      halmd_mdsim_gpu
      halmd_algorithm_gpu
      halmd_utility_gpu
    )
  endif()
  target_link_libraries(test_unit_mdsim_potentials_external_planar_wall
    halmd_mdsim_host_potentials_external_planar_wall
    halmd_mdsim
    ${HALMD_TEST_LIBRARIES}
  )
  add_test(unit/mdsim/potentials/external/planar_wall/host
    test_unit_mdsim_potentials_external_planar_wall --run_test=planar_wall_host --log_level=test_suite
  )
  if(HALMD_WITH_GPU)
    if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
      add_test(unit/mdsim/potentials/external/planar_wall/gpu/float
        test_unit_mdsim_potentials_external_planar_wall --run_test=planar_wall_gpu_float --log_level=test_suite
      )
    endif()
    if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
      add_test(unit/mdsim/potentials/external/planar_wall/gpu/dsfloat
        test_unit_mdsim_potentials_external_planar_wall --run_test=planar_wall_gpu_dsfloat --log_level=test_suite
      )
    endif()
  endif()
endif()
