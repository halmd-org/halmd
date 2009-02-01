if(WITH_CUDA)
  set(ljgpu_VARIANT_CELL_SUMMATION_ORDER TRUE CACHE BOOL
    "Use opposite cell summation order")
  if(ljgpu_VARIANT_CELL_SUMMATION_ORDER)
    add_definitions(-DUSE_CELL_SUMMATION_ORDER)
  endif(ljgpu_VARIANT_CELL_SUMMATION_ORDER)

  set(ljgpu_VARIANT_CELL_DSFUN TRUE CACHE BOOL
    "Use double-single precision functions in cell summation")
  if(ljgpu_VARIANT_CELL_DSFUN)
    add_definitions(-DUSE_CELL_DSFUN)
  endif(ljgpu_VARIANT_CELL_DSFUN)

  set(ljgpu_VARIANT_HILBERT_ORDER TRUE CACHE BOOL
    "Use Hilbert space-filling curve particle ordering")
  if(ljgpu_VARIANT_HILBERT_ORDER)
    add_definitions(-DUSE_HILBERT_ORDER)
  endif(ljgpu_VARIANT_HILBERT_ORDER)

  set(ljgpu_VARIANT_HILBERT_ALT_3D FALSE CACHE BOOL
    "Use alternative 3D Hilbert curve vertex rules")
  if(ljgpu_VARIANT_HILBERT_ALT_3D)
    add_definitions(-DUSE_HILBERT_ALT_3D)
  endif(ljgpu_VARIANT_HILBERT_ALT_3D)
endif(WITH_CUDA)

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
  -Wl,-whole-archive
  util
  sample
  -Wl,-no-whole-archive
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${HDF5_CPP_LIBRARY}
  dl
  )
install(TARGETS ljgpu
  RUNTIME DESTINATION bin
  )

set(ljgpu_impl
    ljfluid_impl_host
    hardsphere_impl
    )
if(WITH_CUDA)
  set(ljgpu_impl
    ${ljgpu_impl}
    ljfluid_impl_gpu_square
    ljfluid_impl_gpu_cell
    ljfluid_impl_gpu_neighbour
    )
endif(WITH_CUDA)

foreach(impl ${ljgpu_impl})
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
    COMPILE_DEFINITIONS "MDSIM_BACKEND=${backend};MDSIM_IMPL=${impl}"
    )
  install(TARGETS ljgpu_${target}
    LIBRARY DESTINATION lib
    )
endforeach(impl)

target_link_libraries(ljgpu_host
  ${GSL_LIBRARIES}
  rt
  )
target_link_libraries(ljgpu_hardsphere
  ${GSL_LIBRARIES}
  rt
  )

if(WITH_CUDA)
  target_link_libraries(ljgpu_gpu_square
    ${CUDA_LIBRARIES}
    gpu_boltzmann
    gpu_lattice
    gpu_ljfluid_square
    gpu_radix_sort
    gpu_rand48
    gpu_reduce
    )
  target_link_libraries(ljgpu_gpu_cell
    ${CUDA_LIBRARIES}
    gpu_boltzmann
    gpu_lattice
    gpu_ljfluid_cell
    gpu_radix_sort
    gpu_rand48
    gpu_reduce
    )
  target_link_libraries(ljgpu_gpu_neighbour
    ${CUDA_LIBRARIES}
    gpu_boltzmann
    gpu_hilbert
    gpu_lattice
    gpu_ljfluid_nbr
    gpu_radix_sort
    gpu_rand48
    gpu_reduce
    )
endif(WITH_CUDA)

add_subdirectory(algorithm)
add_subdirectory(mdsim)
add_subdirectory(rng)
add_subdirectory(sample)
add_subdirectory(test)
add_subdirectory(util)
