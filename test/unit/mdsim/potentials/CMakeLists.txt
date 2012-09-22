add_executable(test_unit_mdsim_potentials_lennard_jones
  lennard_jones.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_potentials_lennard_jones
    halmd_mdsim_gpu_potentials
    halmd_mdsim_gpu
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_potentials_lennard_jones
  halmd_mdsim_host_potentials
  halmd_mdsim
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/potentials/lennard_jones/host
  test_unit_mdsim_potentials_lennard_jones --run_test=lennard_jones_host --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/mdsim/potentials/lennard_jones/gpu
    test_unit_mdsim_potentials_lennard_jones --run_test=lennard_jones_gpu --log_level=test_suite
  )
endif()

add_executable(test_unit_mdsim_potentials_modified_lennard_jones
  modified_lennard_jones.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_potentials_modified_lennard_jones
    halmd_mdsim_gpu_potentials
    halmd_mdsim_gpu
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_potentials_modified_lennard_jones
  halmd_mdsim_host_potentials
  halmd_mdsim
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/potentials/modified_lennard_jones/host
  test_unit_mdsim_potentials_modified_lennard_jones --run_test=modified_lennard_jones_host --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/mdsim/potentials/modified_lennard_jones/gpu
    test_unit_mdsim_potentials_modified_lennard_jones --run_test=modified_lennard_jones_gpu --log_level=test_suite
  )
endif()

add_executable(test_unit_mdsim_potentials_power_law
  power_law.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_potentials_power_law
    halmd_mdsim_gpu_potentials
    halmd_mdsim_gpu
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_potentials_power_law
  halmd_mdsim_host_potentials
  halmd_mdsim
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/potentials/power_law/host
  test_unit_mdsim_potentials_power_law --run_test=power_law_host --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/mdsim/potentials/power_law/gpu
    test_unit_mdsim_potentials_power_law --run_test=power_law_gpu --log_level=test_suite
  )
endif()

add_executable(test_unit_mdsim_potentials_power_law_with_core
  power_law_with_core.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_potentials_power_law_with_core
    halmd_mdsim_gpu_potentials
    halmd_mdsim_gpu
    halmd_algorithm_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_potentials_power_law_with_core
  halmd_mdsim_host_potentials
  halmd_mdsim
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/potentials/power_law_with_core/host
  test_unit_mdsim_potentials_power_law_with_core --run_test=power_law_with_core_host --log_level=test_suite
)
if(HALMD_WITH_GPU)
  add_test(unit/mdsim/potentials/power_law_with_core/gpu
    test_unit_mdsim_potentials_power_law_with_core --run_test=power_law_with_core_gpu --log_level=test_suite
  )
endif()
