if(HALMD_USE_STATIC_LIBS)
  set(HALMD_BACKEND_EXECUTABLES TRUE CACHE BOOL
    "Build separate executable for each backend" FORCE)
else(HALMD_USE_STATIC_LIBS)
  set(HALMD_BACKEND_EXECUTABLES FALSE CACHE BOOL
    "Build separate executable for each backend")
endif(HALMD_USE_STATIC_LIBS)
if(HALMD_BACKEND_EXECUTABLES)
  add_definitions(-DBACKEND_EXECUTABLES)
endif(HALMD_BACKEND_EXECUTABLES)

set(HALMD_VARIANT_HILBERT_ORDER TRUE CACHE BOOL
  "Use Hilbert space-filling curve particle ordering")
if(HALMD_VARIANT_HILBERT_ORDER)
  add_definitions(-DUSE_HILBERT_ORDER)
endif(HALMD_VARIANT_HILBERT_ORDER)

set(HALMD_VARIANT_HILBERT_ALT_3D FALSE CACHE BOOL
  "Use alternative 3D Hilbert curve vertex rules")
if(HALMD_VARIANT_HILBERT_ALT_3D)
  add_definitions(-DUSE_HILBERT_ALT_3D)
endif(HALMD_VARIANT_HILBERT_ALT_3D)

if(WITH_CUDA)
  set(HALMD_VARIANT_CELL_SUMMATION_ORDER TRUE CACHE BOOL
    "Use opposite cell summation order")
  if(HALMD_VARIANT_CELL_SUMMATION_ORDER)
    add_definitions(-DUSE_CELL_SUMMATION_ORDER)
  endif(HALMD_VARIANT_CELL_SUMMATION_ORDER)

  set(HALMD_VARIANT_FORCE_DSFUN TRUE CACHE BOOL
    "Use double-single precision functions in cell summation")
  if(HALMD_VARIANT_FORCE_DSFUN)
    add_definitions(-DUSE_FORCE_DSFUN)
  endif(HALMD_VARIANT_FORCE_DSFUN)

  set(HALMD_VARIANT_VERLET_DSFUN TRUE CACHE BOOL
    "Use double-single precision functions in Verlet integrator")
  if(HALMD_VARIANT_VERLET_DSFUN)
    add_definitions(-DUSE_VERLET_DSFUN)
  endif(HALMD_VARIANT_VERLET_DSFUN)

  set(HALMD_DEVICE_SCALE "3" CACHE STRING
    "Scale/size of the CUDA device (try to reduce in case of insufficient resources)")
  add_definitions(-DDEVICE_SCALE=${HALMD_DEVICE_SCALE})

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
set(HALMD_VARIANT_HOST_SINGLE_PRECISION FALSE CACHE BOOL
  "Use single-precision math in host implementation (requires SSE)")
if(HALMD_VARIANT_HOST_SINGLE_PRECISION)
  add_definitions(-DUSE_HOST_SINGLE_PRECISION)
endif(HALMD_VARIANT_HOST_SINGLE_PRECISION)

# Generate version header file
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/version.h.in
  ${CMAKE_CURRENT_BINARY_DIR}/version.h
  ESCAPE_QUOTES @ONLY
  )

set(HALMD_MODULES
  halmd_io
  halmd_io_profile_writers
  halmd_io_trajectory
  halmd_io_trajectory_readers
  halmd_io_trajectory_writers
  halmd_mdsim
  halmd_mdsim_host
  halmd_mdsim_host_forces
  halmd_mdsim_host_integrator
  halmd_mdsim_host_position
  halmd_mdsim_host_sample
  halmd_mdsim_host_sort
  halmd_mdsim_host_velocities
  halmd_mdsim_samples_host
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
    halmd_random_gpu
    halmd_utility_gpu
    )
endif(WITH_CUDA)

set(HALMD_COMMON_LIBRARIES
  -Wl,-whole-archive
  ${HALMD_MODULES}
  -Wl,-no-whole-archive
  ${Boost_FILESYSTEM_LIBRARY}
  ${Boost_SYSTEM_LIBRARY}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

if(HALMD_USE_STATIC_LIBS)
  set(HALMD_COMMON_LIBRARIES
    ${HALMD_COMMON_LIBRARIES}
    -static-libgcc
    -Wl,-Bstatic
    )
endif(HALMD_USE_STATIC_LIBS)

set(HALMD_COMMON_LIBRARIES
  ${HALMD_COMMON_LIBRARIES}
  rt
  dl
  pthread
  z
  )

if(NOT HALMD_BACKEND_EXECUTABLES)
  add_executable(halmd
    main.cpp
    script.cpp
    )
  target_link_libraries(halmd
    ${HALMD_COMMON_LIBRARIES}
    )
  install(TARGETS halmd
    RUNTIME DESTINATION bin
    )
endif(NOT HALMD_BACKEND_EXECUTABLES)

add_subdirectory(algorithm)
add_subdirectory(io)
add_subdirectory(mdsim)
add_subdirectory(random)
add_subdirectory(test)
add_subdirectory(utility)
