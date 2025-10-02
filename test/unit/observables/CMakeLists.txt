# ssf module
add_executable(test_unit_observables_ssf
  ssf.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_observables_ssf
    halmd_mdsim_gpu_particle_groups
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_observables_ssf
  halmd_mdsim_host_particle_groups
  halmd_mdsim_host_positions
  halmd_mdsim_host
  halmd_mdsim
  halmd_observables_host
  halmd_observables_utility
  halmd_observables
  halmd_random_host
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/observables/ssf/host/2d
  test_unit_observables_ssf --run_test=ssf_host_2d --log_level=test_suite
)
add_test(unit/observables/ssf/host/3d
  test_unit_observables_ssf --run_test=ssf_host_3d --log_level=test_suite
)
if(HALMD_WITH_GPU)
  if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
    halmd_add_gpu_test(unit/observables/ssf/gpu/float/2d
      test_unit_observables_ssf --run_test=ssf_gpu_float_2d --log_level=test_suite
    )
    halmd_add_gpu_test(unit/observables/ssf/gpu/float/3d
      test_unit_observables_ssf --run_test=ssf_gpu_float_3d --log_level=test_suite
    )
  endif()
  if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
    halmd_add_gpu_test(unit/observables/ssf/gpu/dsfloat/2d
      test_unit_observables_ssf --run_test=ssf_gpu_dsfloat_2d --log_level=test_suite
    )
    halmd_add_gpu_test(unit/observables/ssf/gpu/dsfloat/3d
      test_unit_observables_ssf --run_test=ssf_gpu_dsfloat_3d --log_level=test_suite
    )
  endif()
endif()

# phase space sampler
add_executable(test_unit_observables_phase_space
  phase_space.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_observables_phase_space
    halmd_algorithm_gpu
    halmd_mdsim_gpu_particle_groups
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu_velocities
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_observables_phase_space
  halmd_mdsim_host_particle_groups
  halmd_mdsim_host_positions
  halmd_mdsim_host_velocities
  halmd_mdsim_host
  halmd_mdsim
  halmd_observables_host
  halmd_random_host
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/observables/phase_space/host/2d
  test_unit_observables_phase_space --run_test=phase_space_host_2d --log_level=test_suite
)
add_test(unit/observables/phase_space/host/3d
  test_unit_observables_phase_space --run_test=phase_space_host_3d --log_level=test_suite
)
if(HALMD_WITH_GPU)
  if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
    halmd_add_gpu_test(unit/observables/phase_space/gpu/host/float/2d
      test_unit_observables_phase_space --run_test=phase_space_gpu_host_float_2d --log_level=test_suite
    )
    halmd_add_gpu_test(unit/observables/phase_space/gpu/host/float/3d
      test_unit_observables_phase_space --run_test=phase_space_gpu_host_float_3d --log_level=test_suite
    )
    halmd_add_gpu_test(unit/observables/phase_space/gpu/gpu/float/2d
      test_unit_observables_phase_space --run_test=phase_space_gpu_gpu_float_2d --log_level=test_suite
    )
    halmd_add_gpu_test(unit/observables/phase_space/gpu/gpu/float/3d
      test_unit_observables_phase_space --run_test=phase_space_gpu_gpu_float_3d --log_level=test_suite
    )
  endif()
  if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
    halmd_add_gpu_test(unit/observables/phase_space/gpu/host/dsfloat/2d
      test_unit_observables_phase_space --run_test=phase_space_gpu_host_dsfloat_2d --log_level=test_suite
    )
    halmd_add_gpu_test(unit/observables/phase_space/gpu/host/dsfloat/3d
      test_unit_observables_phase_space --run_test=phase_space_gpu_host_dsfloat_3d --log_level=test_suite
    )
    halmd_add_gpu_test(unit/observables/phase_space/gpu/gpu/dsfloat/2d
      test_unit_observables_phase_space --run_test=phase_space_gpu_gpu_dsfloat_2d --log_level=test_suite
    )
    halmd_add_gpu_test(unit/observables/phase_space/gpu/gpu/dsfloat/3d
      test_unit_observables_phase_space --run_test=phase_space_gpu_gpu_dsfloat_3d --log_level=test_suite
    )
  endif()
endif()

# Lua tests for module chemical_potential
configure_file(
  chemical_potential.cmake.in
  test_unit_observables_chemical_potential.cmake
  @ONLY
)

add_test(unit/observables/chemical_potential/lua/host/2d
  ${CMAKE_COMMAND} -DDIMENSION=2 -DUSE_HOST=1 -P test_unit_observables_chemical_potential.cmake
)
add_test(unit/observables/chemical_potential/lua/host/3d
  ${CMAKE_COMMAND} -DDIMENSION=3 -DUSE_HOST=1 -P test_unit_observables_chemical_potential.cmake
)
if(HALMD_WITH_GPU)
  add_test(unit/observables/chemical_potential/lua/gpu/2d
    ${CMAKE_COMMAND} -DDIMENSION=2 -P test_unit_observables_chemical_potential.cmake
  )
  add_test(unit/observables/chemical_potential/lua/gpu/3d
    ${CMAKE_COMMAND} -DDIMENSION=3 -P test_unit_observables_chemical_potential.cmake
  )
endif()

add_subdirectory(utility)
