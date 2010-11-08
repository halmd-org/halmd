set(HALMD_MODULES
  halmd_io
  halmd_io_profile
  halmd_io_profile_writers
  halmd_io_statevars
  halmd_io_statevars_writers
  halmd_io_trajectory
  halmd_io_trajectory_readers
  halmd_io_trajectory_writers
  halmd_io_utility
  halmd_mdsim
  halmd_mdsim_host
  halmd_mdsim_host_forces
  halmd_mdsim_host_integrators
  halmd_mdsim_host_position
  halmd_mdsim_host_sample
  halmd_mdsim_host_sort
  halmd_mdsim_host_velocities
  halmd_mdsim_integrators
  halmd_mdsim_samples_host
  halmd_observables
  halmd_observables_host
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

# If this is a build from the git repository source, we need to generate
# version.h from the current git commit. This should be repeated on every
# invocation of make (that involves dependent source files), so we add a
# phony build rule and let the generator rule depend on it. Both rules
# have empty comments, as make executes them multiple times.
if(GIT_EXECUTABLE)
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/version.h.in_phony
    COMMENT ""
    COMMAND ""
  )
  set_source_files_properties(
    ${CMAKE_CURRENT_BINARY_DIR}/version.h.in_phony PROPERTIES
    SYMBOLIC TRUE
  )
  # The script version.cmake generates version.h with configure_file, which
  # will only regenerate the destination file if the content has changed,
  # so a source rebuild will only be triggered if the git commit changes.
  add_custom_command(
    OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/version.h
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/version.h.in_phony
    COMMENT ""
    COMMAND ${CMAKE_COMMAND}
      -D "INPUT_FILE=${CMAKE_CURRENT_SOURCE_DIR}/version.h.in"
      -D "OUTPUT_FILE=${CMAKE_CURRENT_BINARY_DIR}/version.h"
      -P cmake/version.cmake
  )
  set_source_files_properties(
    ${CMAKE_CURRENT_BINARY_DIR}/version.h PROPERTIES
    GENERATED TRUE
  )
  add_subdirectory(cmake)
endif(GIT_EXECUTABLE)

if(NOT HALMD_BACKEND_EXECUTABLES)
  add_executable(halmd
    main.cpp
    options.cpp
    sampler.cpp
    script.cpp
    version.h
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
