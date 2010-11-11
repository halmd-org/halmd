add_library(halmd_observables_gpu
  thermodynamics.cpp
  thermodynamics_kernel.cu
  trajectory.cpp
  trajectory_kernel.cu
)
add_subdirectory(samples)
