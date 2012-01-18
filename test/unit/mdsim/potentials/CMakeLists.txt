add_executable(test_unit_mdsim_potentials_lennard_jones
  lennard_jones.cpp
)
target_link_libraries(test_unit_mdsim_potentials_lennard_jones
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/potentials/lennard_jones/host
  test_unit_mdsim_potentials_lennard_jones --run_test=lennard_jones_host --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/potentials/lennard_jones/gpu
    test_unit_mdsim_potentials_lennard_jones --run_test=lennard_jones_gpu --log_level=test_suite
  )
endif(WITH_CUDA)

add_executable(test_unit_mdsim_potentials_modified_lennard_jones
  modified_lennard_jones.cpp
)
target_link_libraries(test_unit_mdsim_potentials_modified_lennard_jones
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/potentials/modified_lennard_jones/host
  test_unit_mdsim_potentials_modified_lennard_jones --run_test=modified_lennard_jones_host --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/potentials/modified_lennard_jones/gpu
    test_unit_mdsim_potentials_modified_lennard_jones --run_test=modified_lennard_jones_gpu --log_level=test_suite
  )
endif(WITH_CUDA)

add_executable(test_unit_mdsim_potentials_power_law_with_core
  power_law_with_core.cpp
)
target_link_libraries(test_unit_mdsim_potentials_power_law_with_core
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/potentials/power_law_with_core/host
  test_unit_mdsim_potentials_power_law_with_core --run_test=power_law_with_core_host --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/mdsim/potentials/power_law_with_core/gpu
    test_unit_mdsim_potentials_power_law_with_core --run_test=power_law_with_core_gpu --log_level=test_suite
  )
endif(WITH_CUDA)
