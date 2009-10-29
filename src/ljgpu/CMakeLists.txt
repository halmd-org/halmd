if(ljgpu_USE_STATIC_LIBS)
  set(ljgpu_BACKEND_EXECUTABLES TRUE CACHE BOOL
    "Build separate executable for each backend" FORCE)
else(ljgpu_USE_STATIC_LIBS)
  set(ljgpu_BACKEND_EXECUTABLES FALSE CACHE BOOL
    "Build separate executable for each backend")
endif(ljgpu_USE_STATIC_LIBS)
if(ljgpu_BACKEND_EXECUTABLES)
  add_definitions(-DBACKEND_EXECUTABLES)
endif(ljgpu_BACKEND_EXECUTABLES)

if(WITH_CUDA)
  set(ljgpu_VARIANT_CELL_SUMMATION_ORDER TRUE CACHE BOOL
    "Use opposite cell summation order")
  if(ljgpu_VARIANT_CELL_SUMMATION_ORDER)
    add_definitions(-DUSE_CELL_SUMMATION_ORDER)
  endif(ljgpu_VARIANT_CELL_SUMMATION_ORDER)

  set(ljgpu_VARIANT_FORCE_DSFUN TRUE CACHE BOOL
    "Use double-single precision functions in cell summation")
  if(ljgpu_VARIANT_FORCE_DSFUN)
    add_definitions(-DUSE_FORCE_DSFUN)
  endif(ljgpu_VARIANT_FORCE_DSFUN)

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

  set(ljgpu_VARIANT_VERLET_DSFUN TRUE CACHE BOOL
    "Use double-single precision functions in Verlet integrator")
  if(ljgpu_VARIANT_VERLET_DSFUN)
    add_definitions(-DUSE_VERLET_DSFUN)
  endif(ljgpu_VARIANT_VERLET_DSFUN)

  set(ljgpu_DEVICE_SCALE "3" CACHE STRING
    "Scale/size of the CUDA device (try to reduce in case of insufficient resources)")
  add_definitions(-DDEVICE_SCALE=${ljgpu_DEVICE_SCALE})

endif(WITH_CUDA)

#
# The following option only works on x86-64 by default, which always uses the
# SSE instruction set for floating-point math. On i386, the x87 floating-point
# unit provides 80-bit extended double precision math internally, causing
# excess precision even if values are stored in single-precision.
#
# See the -mfpmath option in the gcc manpage for details, and
#
# Deterministic cross-platform floating point arithmetics
# http://www.christian-seiler.de/projekte/fpmath/
#
set(ljgpu_VARIANT_HOST_SINGLE_PRECISION FALSE CACHE BOOL
  "Use single-precision math in host implementation (requires SSE)")
if(ljgpu_VARIANT_HOST_SINGLE_PRECISION)
  add_definitions(-DUSE_HOST_SINGLE_PRECISION)
endif(ljgpu_VARIANT_HOST_SINGLE_PRECISION)

# Generate version header file
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.h
  ESCAPE_QUOTES @ONLY
  )

set(ljgpu_COMMON_LIBRARIES
  -Wl,-whole-archive
  util
  sample
  -Wl,-no-whole-archive
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${HDF5_CPP_LIBRARY}
  ${HDF5_LIBRARY}
  )

if(ljgpu_USE_STATIC_LIBS)
  set(ljgpu_COMMON_LIBRARIES
    ${ljgpu_COMMON_LIBRARIES}
    -static-libgcc
    -Wl,-Bstatic
    )
endif(ljgpu_USE_STATIC_LIBS)

set(ljgpu_COMMON_LIBRARIES
  ${ljgpu_COMMON_LIBRARIES}
  rt
  dl
  pthread
  z
  )

if(NOT ljgpu_BACKEND_EXECUTABLES)
  add_executable(ljgpu
    main.cpp
    options.cpp
    )
  target_link_libraries(ljgpu
    ${ljgpu_COMMON_LIBRARIES}
    )
  install(TARGETS ljgpu
    RUNTIME DESTINATION bin
    )
endif(NOT ljgpu_BACKEND_EXECUTABLES)

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
  string(REGEX REPLACE "_impl[_a-z]*$" "" class ${impl})
  string(REGEX REPLACE "(^ljfluid_impl_|_impl)" "" backend ${impl})

  set(ljgpu_COMPILE_DEFINITIONS
    MDSIM_CLASS=${class}
    MDSIM_IMPL=${impl}
    MDSIM_BACKEND="${backend}"
    )

  if(ljgpu_BACKEND_EXECUTABLES)
    add_executable(ljgpu_${backend}
      main.cpp
      options.cpp
      mdlib.cpp
      )
    target_link_libraries(ljgpu_${backend}
      ${ljgpu_COMMON_LIBRARIES}
      )
    install(TARGETS ljgpu_${backend}
      RUNTIME DESTINATION bin
      )
  else(ljgpu_BACKEND_EXECUTABLES)
    add_library(ljgpu_${backend} MODULE
      mdlib.cpp
      )
    install(TARGETS ljgpu_${backend}
      LIBRARY DESTINATION lib
      )
  endif(ljgpu_BACKEND_EXECUTABLES)

  set_target_properties(ljgpu_${backend} PROPERTIES
    COMPILE_DEFINITIONS "${ljgpu_COMPILE_DEFINITIONS}"
    )
  if(ljgpu_USE_STATIC_LIBS)
    set_target_properties(ljgpu_${backend} PROPERTIES
      LINK_SEARCH_END_STATIC TRUE
      )
  endif(ljgpu_USE_STATIC_LIBS)
endforeach(impl)

target_link_libraries(ljgpu_host
  ${GSL_LIBRARIES}
  )
target_link_libraries(ljgpu_hardsphere
  ${GSL_LIBRARIES}
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
    gpu_virial
    )
  target_link_libraries(ljgpu_gpu_cell
    ${CUDA_LIBRARIES}
    gpu_boltzmann
    gpu_lattice
    gpu_ljfluid_cell
    gpu_radix_sort
    gpu_rand48
    gpu_reduce
    gpu_virial
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
    gpu_virial
    )
endif(WITH_CUDA)

add_subdirectory(algorithm)
add_subdirectory(mdsim)
add_subdirectory(rng)
add_subdirectory(sample)
add_subdirectory(test)
add_subdirectory(util)
