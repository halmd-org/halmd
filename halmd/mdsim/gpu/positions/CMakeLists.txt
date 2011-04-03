add_library(halmd_mdsim_gpu_positions_lattice
  lattice.cpp
  lattice_kernel.cu
)
add_library(halmd_mdsim_gpu_positions_phase_space
  phase_space.cpp
  phase_space_kernel.cu
)
