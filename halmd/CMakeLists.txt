# Generate version header file
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.h
  ESCAPE_QUOTES @ONLY
  )

set(HALMD_MODULES
  halmd_io
  halmd_io_profile
  halmd_io_profile_writers
  halmd_io_statevars
  halmd_io_statevars_writers
  halmd_io_trajectory
  halmd_io_trajectory_readers
  halmd_io_trajectory_writers
  halmd_mdsim
  halmd_mdsim_host
  halmd_mdsim_host_forces
  halmd_mdsim_host_integrators
  halmd_mdsim_host_position
  halmd_mdsim_host_sample
  halmd_mdsim_host_sort
  halmd_mdsim_host_velocities
  halmd_mdsim_samples_host
  halmd_observables
  halmd_observables_host
  halmd_random
  halmd_random_host
  halmd_utility
  )
if(WITH_CUDA)
  set(HALMD_MODULES
    ${HALMD_MODULES}
    halmd_mdsim_gpu
    halmd_mdsim_gpu_forces
    halmd_mdsim_gpu_integrators
    halmd_mdsim_gpu_position
    halmd_mdsim_gpu_sample
    halmd_mdsim_gpu_sort
    halmd_mdsim_gpu_velocities
    halmd_mdsim_samples_gpu
    halmd_observables_gpu
    halmd_random_gpu
    halmd_utility_gpu
    )
endif(WITH_CUDA)

if(NOT HALMD_BACKEND_EXECUTABLES)
  add_executable(halmd
    main.cpp
    options.cpp
    sampler.cpp
    script.cpp
    )
  target_link_libraries(halmd
    -Wl,-whole-archive
    ${HALMD_MODULES}
    -Wl,-no-whole-archive
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
