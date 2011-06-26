# thermodynamics module
add_executable(test_unit_observables_thermodynamics
  thermodynamics.cpp
)
target_link_libraries(test_unit_observables_thermodynamics
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/observables/thermodynamics/host/2d
  test_unit_observables_thermodynamics --run_test=lennard_jones_fluid_host_2d --log_level=test_suite
)
add_test(unit/observables/thermodynamics/host/3d
  test_unit_observables_thermodynamics --run_test=lennard_jones_fluid_host_3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/observables/thermodynamics/gpu/2d
    test_unit_observables_thermodynamics --run_test=lennard_jones_fluid_gpu_2d --log_level=test_suite
  )
  add_test(unit/observables/thermodynamics/gpu/3d
    test_unit_observables_thermodynamics --run_test=lennard_jones_fluid_gpu_3d --log_level=test_suite
  )
endif(WITH_CUDA)

# ssf module
add_executable(test_unit_observables_ssf
  ssf.cpp
)
target_link_libraries(test_unit_observables_ssf
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/observables/ssf/host/2d
  test_unit_observables_ssf --run_test=ssf_host_2d --log_level=test_suite
)
add_test(unit/observables/ssf/host/3d
  test_unit_observables_ssf --run_test=ssf_host_3d --log_level=test_suite
)
if(WITH_CUDA)
  add_test(unit/observables/ssf/gpu/2d
    test_unit_observables_ssf --run_test=ssf_gpu_2d --log_level=test_suite
  )
  add_test(unit/observables/ssf/gpu/3d
    test_unit_observables_ssf --run_test=ssf_gpu_3d --log_level=test_suite
  )
endif(WITH_CUDA)
