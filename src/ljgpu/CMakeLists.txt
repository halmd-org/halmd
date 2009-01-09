add_subdirectory(algorithm)
add_subdirectory(ljfluid)
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
  mdsample
  mdutil
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${HDF5_CPP_LIBRARY}
  ${CUDA_LIBRARIES}
  gpu_hilbert
  gpu_lattice
  gpu_ljfluid_cell
  gpu_ljfluid_nbr
  gpu_ljfluid_square
  gpu_radix_sort
  gpu_rand48
  gpu_reduce
  ${GSL_LIBRARIES}
  rt
  )
install(TARGETS ljgpu
  RUNTIME DESTINATION bin
  )
