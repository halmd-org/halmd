##
# HALMD C++ modules
#
# Base classes must be ordered *before* derived classes, otherwise
# Luabind will throw an assertion error or cause a segmentation fault.
#
set(HALMD_MODULES
  halmd_io_logger
  halmd_io_profiling_writer
  halmd_io_profiling_writers_h5md
  halmd_io_profiling_writers_log
  halmd_io_statevars_writer
  halmd_io_statevars_writers_hdf5
  halmd_io_trajectory_reader
  halmd_io_trajectory_readers_h5md
  halmd_io_trajectory_readers_halmd_0_1_x
  halmd_io_trajectory_writer
  halmd_io_trajectory_writers_h5md
  halmd_mdsim_box
  halmd_mdsim_core
  halmd_mdsim_force
  halmd_mdsim_integrator
  halmd_mdsim_integrators_nvt
  halmd_mdsim_neighbour
  halmd_mdsim_particle
  halmd_mdsim_position
  halmd_mdsim_sort
  halmd_mdsim_velocity
  halmd_mdsim_host_force
  halmd_mdsim_host_forces_lennard_jones
  halmd_mdsim_host_forces_morse
  halmd_mdsim_host_forces_power_law
  halmd_mdsim_host_forces_smooth
  halmd_mdsim_host_forces_zero
  halmd_mdsim_host_integrators_verlet
  halmd_mdsim_host_integrators_verlet_nvt_andersen
  halmd_mdsim_host_neighbour
  halmd_mdsim_host_particle
  halmd_mdsim_host_positions_lattice
  halmd_mdsim_host_positions_phase_space
  halmd_mdsim_host_sorts_hilbert
  halmd_mdsim_host_velocity
  halmd_mdsim_host_velocities_boltzmann
  halmd_mdsim_host_velocities_phase_space
  halmd_runner
  halmd_observables_density_mode
  halmd_observables_phase_space
  halmd_observables_sampler
  halmd_observables_samples_density_mode
  halmd_observables_ssf
  halmd_observables_thermodynamics
  halmd_observables_utility_wavevector
  halmd_observables_host_density_mode
  halmd_observables_host_mean_quartic_displacement
  halmd_observables_host_mean_square_displacement
  halmd_observables_host_phase_space
  halmd_observables_host_samples_phase_space
  halmd_observables_host_thermodynamics
  halmd_observables_host_velocity_autocorrelation
  halmd_random_random
  halmd_random_host_random
  halmd_utility_lua_any_converter
  halmd_utility_lua_hdf5
  halmd_utility_lua_program_options
  halmd_utility_lua_signal
  halmd_utility_lua_ublas
  halmd_utility_options_parser
  halmd_utility_profiler
)
if(WITH_CUDA)
  set(HALMD_MODULES
    ${HALMD_MODULES}
    halmd_mdsim_gpu_force
    halmd_mdsim_gpu_forces_lennard_jones
    halmd_mdsim_gpu_forces_morse
    halmd_mdsim_gpu_forces_zero
    halmd_mdsim_gpu_integrators_verlet
    halmd_mdsim_gpu_integrators_verlet_nvt_andersen
    halmd_mdsim_gpu_neighbour
    halmd_mdsim_gpu_particle
    halmd_mdsim_gpu_positions_lattice
    halmd_mdsim_gpu_positions_phase_space
    halmd_mdsim_gpu_sorts_hilbert
    halmd_mdsim_gpu_velocity
    halmd_mdsim_gpu_velocities_boltzmann
    halmd_mdsim_gpu_velocities_phase_space
    halmd_observables_gpu_density_mode
    halmd_observables_gpu_phase_space
    halmd_observables_gpu_samples_phase_space
    halmd_observables_gpu_thermodynamics
    halmd_random_gpu_random
    halmd_utility_gpu_device
  )
endif(WITH_CUDA)

add_library(halmd_runner
  runner.cpp
)

add_subdirectory(algorithm)
add_subdirectory(io)
add_subdirectory(mdsim)
add_subdirectory(observables)
add_subdirectory(random)
add_subdirectory(utility)

##
# Generate Lua binding function calls from module target sources
#
set(modules_hpp)
set(modules_cpp)
foreach(target ${HALMD_MODULES})
  get_target_property(sources ${target} SOURCES)
  get_target_property(location ${target} LOCATION)
  get_filename_component(binary_dir "${location}" PATH)
  file(RELATIVE_PATH subdir "${CMAKE_CURRENT_BINARY_DIR}" "${binary_dir}")
  set(source_dir "${CMAKE_CURRENT_SOURCE_DIR}/${subdir}")
  foreach(source ${sources})
    if(NOT IS_ABSOLUTE "${source}")
      set(source "${source_dir}/${source}")
    endif()
    file(STRINGS "${source}" line REGEX "HALMD_LUA_API")
    if(line MATCHES "^HALMD_LUA_API int ([a-z0-9_]+)\\(lua_State\\* L\\)$")
      set(modules_hpp "${modules_hpp}\n${line};")
      set(modules_cpp "${modules_cpp}\n    ${CMAKE_MATCH_1}(L);")
    elseif(line)
      message(SEND_ERROR "malformed HALMD_LUA_API function definition in ${source}")
    endif()
  endforeach()
endforeach()
string(STRIP "${modules_hpp}" HALMD_MODULES_LUAOPEN_DECLARE)
string(STRIP "${modules_cpp}" HALMD_MODULES_LUAOPEN_CALL)

configure_file(modules.hpp.in modules.hpp @ONLY)
configure_file(modules.cpp.in modules.cpp @ONLY)
configure_file(version.h.in version.h ESCAPE_QUOTES @ONLY)

add_executable(halmd
  main.cpp
  "${CMAKE_CURRENT_BINARY_DIR}/modules.cpp"
  script.cpp
)
target_link_libraries(halmd
  ${HALMD_MODULES}
  halmd_utility_hostname
  ${HALMD_COMMON_LIBRARIES}
)
install(TARGETS halmd
  RUNTIME DESTINATION bin
)
