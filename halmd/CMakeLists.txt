set(HALMD_MODULES
  halmd_io
  halmd_io_profiling
  halmd_io_profiling_writers
  halmd_io_statevars
  halmd_io_statevars_writers
  halmd_io_trajectory
  halmd_io_trajectory_readers
  halmd_io_trajectory_writers
  halmd_mdsim
  halmd_mdsim_host
  halmd_mdsim_host_forces
  halmd_mdsim_host_integrators
  halmd_mdsim_host_positions
  halmd_mdsim_host_sorts
  halmd_mdsim_host_velocities
  halmd_mdsim_integrators
  halmd_observables
  halmd_observables_host
  halmd_observables_host_samples
  halmd_observables_samples
  halmd_observables_utility
  halmd_random
  halmd_random_host
  halmd_utility
  halmd_utility_lua_wrapper
  )
if(WITH_CUDA)
  set(HALMD_MODULES
    ${HALMD_MODULES}
    halmd_mdsim_gpu
    halmd_mdsim_gpu_forces
    halmd_mdsim_gpu_integrators
    halmd_mdsim_gpu_positions
    halmd_mdsim_gpu_sorts
    halmd_mdsim_gpu_velocities
    halmd_observables_gpu
    halmd_observables_gpu_samples
    halmd_random_gpu
    halmd_utility_gpu
    )
endif(WITH_CUDA)

configure_file(version.h.in version.h ESCAPE_QUOTES @ONLY)

if(NOT HALMD_BACKEND_EXECUTABLES)
  add_executable(halmd
    main.cpp
    modules.cpp
    runner.cpp
    script.cpp
  )
  target_link_libraries(halmd
    -Wl,-whole-archive
    ${HALMD_MODULES}
    -Wl,-no-whole-archive
    halmd_utility_program_options
    ${HALMD_COMMON_LIBRARIES}
    )
  install(TARGETS halmd
    RUNTIME DESTINATION bin
    )
endif(NOT HALMD_BACKEND_EXECUTABLES)

add_subdirectory(algorithm)
add_subdirectory(io)
add_subdirectory(mdsim)
add_subdirectory(observables)
add_subdirectory(random)
add_subdirectory(utility)
