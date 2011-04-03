add_library(halmd_mdsim_gpu_positions
  lattice.cpp
  lattice_kernel.cu
  phase_space.cpp
  phase_space_kernel.cu
)
