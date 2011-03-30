##
# HALMD C++ modules
#
# Base classes must be ordered *before* derived classes, otherwise
# Luabind will throw an assertion error or cause a segmentation fault.
#
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
  halmd_mdsim_integrators
  halmd_mdsim_host
  halmd_mdsim_host_forces
  halmd_mdsim_host_integrators
  halmd_mdsim_host_positions
  halmd_mdsim_host_sorts
  halmd_mdsim_host_velocities
  halmd_runner
  halmd_observables
  halmd_observables_host
  halmd_observables_host_samples
  halmd_observables_samples
  halmd_observables_utility
  halmd_random
  halmd_random_host
  halmd_utility
  halmd_utility_lua_wrapper
  halmd_utility_program_options
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
  ${HALMD_COMMON_LIBRARIES}
)
install(TARGETS halmd
  RUNTIME DESTINATION bin
)
