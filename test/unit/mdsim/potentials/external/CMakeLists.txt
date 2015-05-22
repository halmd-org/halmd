if(${HALMD_WITH_external_slit})
  add_executable(test_unit_mdsim_potentials_external_slit
    slit.cpp
  )
  if(HALMD_WITH_GPU)
    target_link_libraries(test_unit_mdsim_potentials_external_slit
      halmd_mdsim_gpu_potentials_external_slit
      halmd_mdsim_gpu
      halmd_algorithm_gpu
      halmd_utility_gpu
    )
  endif()
  target_link_libraries(test_unit_mdsim_potentials_external_slit
    halmd_mdsim_host_potentials_external_slit
    halmd_mdsim
    ${HALMD_TEST_LIBRARIES}
  )
  add_test(unit/mdsim/potentials/external/slit/host
    test_unit_mdsim_potentials_external_slit --run_test=slit_host --log_level=test_suite
  )
  if(HALMD_WITH_GPU)
    add_test(unit/mdsim/potentials/external/slit/gpu
      test_unit_mdsim_potentials_external_slit --run_test=slit_gpu --log_level=test_suite
    )
  endif()
endif()
