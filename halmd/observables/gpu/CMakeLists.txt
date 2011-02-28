add_library(halmd_observables_gpu
  density_modes.cpp
#   density_modes_kernel.cu
  thermodynamics.cpp
  thermodynamics_kernel.cu
  trajectory.cpp
  trajectory_kernel.cu
)
add_subdirectory(samples)
