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

