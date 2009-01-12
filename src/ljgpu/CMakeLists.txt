add_subdirectory(algorithm)
add_subdirectory(mdsim)
add_subdirectory(rng)
add_subdirectory(sample)
add_subdirectory(util)

# Generate version header file
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.h
  ESCAPE_QUOTES @ONLY
  )

add_executable(ljgpu
  main.cpp
  options.cpp
  )
target_link_libraries(ljgpu
  util_h5xx
  util_log
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${HDF5_CPP_LIBRARY}
  ${CUDA_LIBRARIES}
  )
install(TARGETS ljgpu
  RUNTIME DESTINATION bin
  )

foreach(impl
    ljfluid_impl_host
    ljfluid_impl_gpu_square
    ljfluid_impl_gpu_cell
    ljfluid_impl_gpu_neighbour
    hardsphere_impl
    )
  string(REGEX REPLACE "_impl[_a-z]*$" "" backend ${impl})
  if(impl MATCHES "_impl$")
    string(REGEX REPLACE "_impl$" "" target ${impl})
  else(impl MATCHES "_impl$")
    string(REGEX REPLACE "^[a-z]+_impl_" "" target ${impl})
  endif(impl MATCHES "_impl$")

  add_library(ljgpu_${target} MODULE
    mdlib.cpp
    )
  set_target_properties(ljgpu_${target} PROPERTIES
    COMPILE_DEFINITIONS "MDSIM_BACKEND=${backend};MDSIM_IMPL=${impl}")
  target_link_libraries(ljgpu_${target}
    sample
    util_h5param
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${HDF5_CPP_LIBRARY}
    )
  install(TARGETS ljgpu_${target}
    LIBRARY DESTINATION bin
    )
endforeach(impl)

target_link_libraries(ljgpu_host
  ${GSL_LIBRARIES}
  rt
  )

target_link_libraries(ljgpu_gpu_square
  ${CUDA_LIBRARIES}
  gpu_lattice
  gpu_ljfluid_square
  gpu_rand48
  gpu_reduce
  )

target_link_libraries(ljgpu_gpu_cell
  ${CUDA_LIBRARIES}
  gpu_lattice
  gpu_ljfluid_cell
  gpu_rand48
  )

target_link_libraries(ljgpu_gpu_neighbour
  ${CUDA_LIBRARIES}
  gpu_hilbert
  gpu_lattice
  gpu_ljfluid_nbr
  gpu_radix_sort
  gpu_rand48
  gpu_reduce
  )

target_link_libraries(ljgpu_hardsphere
  ${GSL_LIBRARIES}
  rt
  )
