# module brownian
add_executable(test_unit_mdsim_integrators_brownian
  brownian.cpp
  apply_negate.cu
)
if(HALMD_WITH_GPU)
    target_link_libraries(test_unit_mdsim_integrators_brownian
   halmd_mdsim_gpu_integrators
   halmd_mdsim_gpu_particle_groups
   halmd_mdsim_gpu_positions
   halmd_mdsim_gpu_orientations
   halmd_mdsim_gpu_velocities
   halmd_mdsim_gpu
   halmd_observables_gpu
   halmd_observables_gpu_dynamics
   halmd_random_gpu
   halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_integrators_brownian
  halmd_mdsim_host_integrators
  halmd_mdsim_host_particle_groups
  halmd_mdsim_host_positions
  halmd_mdsim_host_velocities
  halmd_mdsim_host
  halmd_mdsim
  halmd_observables_dynamics
  halmd_observables_host
  halmd_observables_host_dynamics
  halmd_random_host
  ${HALMD_TEST_LIBRARIES}
  )
  #add_test(unit/mdsim/integrators/brownian/host/3d/overdamped
  #test_unit_mdsim_integrators_brownian --run_test=brownian_host_3d_overdamped --log_level=test_suite
  #)
if(HALMD_WITH_GPU)
  # add_test(unit/mdsim/integrators/brownian/gpu/2d/linear
  # test_unit_mdsim_integrators_brownian --run_test=brownian_gpu_2d_linear --log_level=test_suite
  #)
  #add_test(unit/mdsim/integrators/brownian/gpu/2d/overdamped
  # test_unit_mdsim_integrators_brownian --run_test=brownian_gpu_2d_overdamped --log_level=test_suite
  #)
  add_test(unit/mdsim/integrators/brownian/gpu/3d/overdamped
    test_unit_mdsim_integrators_brownian --run_test=brownian_gpu_3d_overdamped --log_level=test_suite
  )
  #add_test(unit/mdsim/integrators/brownian/gpu/3d/overdamped
  #  test_unit_mdsim_integrators_brownian --run_test=brownian_gpu_3d_overdamped --log_level=test_suite
  #)
endif()

# module euler
cuda_add_executable(test_unit_mdsim_integrators_euler
  euler.cpp
  apply_negate.cu
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_integrators_euler
    halmd_mdsim_gpu_integrators
    halmd_mdsim_gpu_particle_groups
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu_velocities
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_integrators_euler
  halmd_mdsim_host_integrators
  halmd_mdsim_host_particle_groups
  halmd_mdsim_host_positions
  halmd_mdsim_host_velocities
  halmd_mdsim_host
  halmd_mdsim
  halmd_observables_host
  halmd_random_host
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
if(HALMD_WITH_GPU)
  if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/euler/gpu/float/2d/linear
      test_unit_mdsim_integrators_euler --run_test=euler_gpu_float_2d_linear --log_level=test_suite
    )
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/euler/gpu/float/2d/overdamped
      test_unit_mdsim_integrators_euler --run_test=euler_gpu_float_2d_overdamped --log_level=test_suite
    )
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/euler/gpu/float/3d/linear
      test_unit_mdsim_integrators_euler --run_test=euler_gpu_float_3d_linear --log_level=test_suite
    )
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/euler/gpu/float/3d/overdamped
      test_unit_mdsim_integrators_euler --run_test=euler_gpu_float_3d_overdamped --log_level=test_suite
    )
  endif()
  if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/euler/gpu/dsfloat/2d/linear
      test_unit_mdsim_integrators_euler --run_test=euler_gpu_dsfloat_2d_linear --log_level=test_suite
    )
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/euler/gpu/dsfloat/2d/overdamped
      test_unit_mdsim_integrators_euler --run_test=euler_gpu_dsfloat_2d_overdamped --log_level=test_suite
    )
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/euler/gpu/dsfloat/3d/linear
      test_unit_mdsim_integrators_euler --run_test=euler_gpu_dsfloat_3d_linear --log_level=test_suite
    )
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/euler/gpu/dsfloat/3d/overdamped
      test_unit_mdsim_integrators_euler --run_test=euler_gpu_dsfloat_3d_overdamped --log_level=test_suite
    )
  endif()
endif()

# module verlet
add_executable(test_unit_mdsim_integrators_verlet
  verlet.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_integrators_verlet
    halmd_mdsim_gpu_integrators
    halmd_mdsim_gpu_particle_groups
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu_velocities
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_integrators_verlet
  halmd_mdsim_host_particle_groups
  halmd_mdsim_host_integrators
  halmd_mdsim_host_positions
  halmd_mdsim_host_velocities
  halmd_mdsim_host
  halmd_mdsim
  halmd_observables_host
  halmd_observables
  halmd_random_host
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/integrators/verlet/host/2d
  test_unit_mdsim_integrators_verlet --run_test=ideal_gas_host_2d --log_level=test_suite
)
add_test(unit/mdsim/integrators/verlet/host/3d
  test_unit_mdsim_integrators_verlet --run_test=ideal_gas_host_3d --log_level=test_suite
)
if(HALMD_WITH_GPU)
  if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet/gpu/float/2d
      test_unit_mdsim_integrators_verlet --run_test=ideal_gas_gpu_float_2d --log_level=test_suite
    )
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet/gpu/float/3d
      test_unit_mdsim_integrators_verlet --run_test=ideal_gas_gpu_float_3d --log_level=test_suite
    )
  endif()
  if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet/gpu/dsfloat/2d
      test_unit_mdsim_integrators_verlet --run_test=ideal_gas_gpu_dsfloat_2d --log_level=test_suite
    )
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet/gpu/dsfloat/3d
      test_unit_mdsim_integrators_verlet --run_test=ideal_gas_gpu_dsfloat_3d --log_level=test_suite
    )
  endif()
endif()

# module verlet_nvt_andersen
add_executable(test_unit_mdsim_integrators_verlet_nvt_andersen
  verlet_nvt_andersen.cpp
)
if(HALMD_WITH_GPU)
  target_link_libraries(test_unit_mdsim_integrators_verlet_nvt_andersen
    halmd_mdsim_gpu_integrators
    halmd_mdsim_gpu_particle_groups
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu_velocities
    halmd_mdsim_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
  )
endif()
target_link_libraries(test_unit_mdsim_integrators_verlet_nvt_andersen
  halmd_mdsim_host_integrators
  halmd_mdsim_host_particle_groups
  halmd_mdsim_host_positions
  halmd_mdsim_host_velocities
  halmd_mdsim_host
  halmd_mdsim
  halmd_observables_host
  halmd_observables
  halmd_random_host
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/integrators/verlet_nvt_andersen/host/2d
  test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen_host_2d --log_level=test_suite
)
add_test(unit/mdsim/integrators/verlet_nvt_andersen/host/3d
  test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen_host_3d --log_level=test_suite
)
if(HALMD_WITH_GPU)
  if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet_nvt_andersen/gpu/float/2d
      test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen_gpu_float_2d --log_level=test_suite
    )
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet_nvt_andersen/gpu/float/3d
      test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen_gpu_float_3d --log_level=test_suite
    )
  endif()
  if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet_nvt_andersen/gpu/dsfloat/2d
      test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen_gpu_dsfloat_2d --log_level=test_suite
    )
    halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet_nvt_andersen/gpu/dsfloat/3d
      test_unit_mdsim_integrators_verlet_nvt_andersen --run_test=verlet_nvt_andersen_gpu_dsfloat_3d --log_level=test_suite
    )
  endif()
endif()

# module verlet_nvt_hoover
if(HALMD_WITH_pair_lennard_jones)
  add_executable(test_unit_mdsim_integrators_verlet_nvt_hoover
    verlet_nvt_hoover.cpp
  )
  if(HALMD_WITH_GPU)
    target_link_libraries(test_unit_mdsim_integrators_verlet_nvt_hoover
      halmd_mdsim_gpu_particle_groups
      halmd_mdsim_gpu_integrators
      halmd_mdsim_gpu_neighbours
      halmd_mdsim_gpu_positions
      halmd_mdsim_gpu_potentials_pair_lennard_jones
      halmd_mdsim_gpu_velocities
      halmd_mdsim_gpu
      halmd_observables_gpu
      halmd_random_gpu
      halmd_utility_gpu
    )
  endif()
  target_link_libraries(test_unit_mdsim_integrators_verlet_nvt_hoover
    halmd_mdsim_host_integrators
    halmd_mdsim_host_neighbours
    halmd_mdsim_host_particle_groups
    halmd_mdsim_host_positions
    halmd_mdsim_host_potentials_pair_lennard_jones
    halmd_mdsim_host_velocities
    halmd_mdsim_host
    halmd_mdsim
    halmd_observables_host
    halmd_observables
    halmd_random_host
    ${HALMD_TEST_LIBRARIES}
  )
  add_test(unit/mdsim/integrators/verlet_nvt_hoover/host/2d
    test_unit_mdsim_integrators_verlet_nvt_hoover --run_test=verlet_nvt_hoover_host_2d --log_level=test_suite
  )
  add_test(unit/mdsim/integrators/verlet_nvt_hoover/host/3d
    test_unit_mdsim_integrators_verlet_nvt_hoover --run_test=verlet_nvt_hoover_host_3d --log_level=test_suite
  )
  if(HALMD_WITH_GPU)
    if(HALMD_VARIANT_GPU_SINGLE_PRECISION)
      halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet_nvt_hoover/gpu/float/2d
        test_unit_mdsim_integrators_verlet_nvt_hoover --run_test=verlet_nvt_hoover_gpu_float_2d --log_level=test_suite
      )
      halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet_nvt_hoover/gpu/float/3d
        test_unit_mdsim_integrators_verlet_nvt_hoover --run_test=verlet_nvt_hoover_gpu_float_3d --log_level=test_suite
      )
    endif()
    if(HALMD_VARIANT_GPU_DOUBLE_SINGLE_PRECISION)
      halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet_nvt_hoover/gpu/dsfloat/2d
        test_unit_mdsim_integrators_verlet_nvt_hoover --run_test=verlet_nvt_hoover_gpu_dsfloat_2d --log_level=test_suite
      )
      halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet_nvt_hoover/gpu/dsfloat/3d
        test_unit_mdsim_integrators_verlet_nvt_hoover --run_test=verlet_nvt_hoover_gpu_dsfloat_3d --log_level=test_suite
      )
    endif()
  endif()
endif()

if(HALMD_VARIANT_HOST_SINGLE_PRECISION)
  set(PARAMETER_TOLERANCE "1e-7")
else()
  set(PARAMETER_TOLERANCE "1e-15")
endif()
# Lua tests for module verlet_nvt_hoover
configure_file(
  verlet_nvt_hoover.cmake.in
  test_unit_mdsim_integrators_verlet_nvt_hoover.cmake
  @ONLY
)

add_test(unit/mdsim/integrators/verlet_nvt_hoover/lua/host/2d
  ${CMAKE_COMMAND} -DDIMENSION=2 -DUSE_HOST=1 -P test_unit_mdsim_integrators_verlet_nvt_hoover.cmake
)
add_test(unit/mdsim/integrators/verlet_nvt_hoover/lua/host/3d
  ${CMAKE_COMMAND} -DDIMENSION=3 -DUSE_HOST=1 -P test_unit_mdsim_integrators_verlet_nvt_hoover.cmake
)
if(HALMD_WITH_GPU)
  halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet_nvt_hoover/lua/gpu/2d
    ${CMAKE_COMMAND} -DDIMENSION=2 -P test_unit_mdsim_integrators_verlet_nvt_hoover.cmake
  )
  halmd_add_gpu_test(NO_MEMCHECK unit/mdsim/integrators/verlet_nvt_hoover/lua/gpu/3d
    ${CMAKE_COMMAND} -DDIMENSION=3 -P test_unit_mdsim_integrators_verlet_nvt_hoover.cmake
  )
endif()
