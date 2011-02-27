add_library(halmd_mdsim_gpu_positions
  lattice.cpp
  lattice_kernel.cu
  trajectory.cpp
  trajectory_kernel.cu
)
