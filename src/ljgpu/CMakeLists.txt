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
  mdutil
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${HDF5_CPP_LIBRARY}
  ${CUDA_LIBRARIES}
  )
install(TARGETS ljgpu
  RUNTIME DESTINATION bin
  )

foreach(LJFLUID_IMPL host gpu_square gpu_cell gpu_neighbour)
  add_library(ljgpu_${LJFLUID_IMPL} MODULE
    mdlib.cpp
    )
  set_target_properties(ljgpu_${LJFLUID_IMPL} PROPERTIES
    COMPILE_DEFINITIONS "LJFLUID_IMPL=ljfluid_impl_${LJFLUID_IMPL}")
  target_link_libraries(ljgpu_${LJFLUID_IMPL}
    mdsample
    mdutil
    ${Boost_FILESYSTEM_LIBRARY}
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${HDF5_CPP_LIBRARY}
    )
  install(TARGETS ljgpu_${LJFLUID_IMPL}
    LIBRARY DESTINATION bin
    )
endforeach(LJFLUID_IMPL)

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
