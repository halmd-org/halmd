halmd_add_library(halmd_mdsim_gpu_positions
  lattice.cpp
  lattice_kernel.cu
)
halmd_add_modules(
  libhalmd_mdsim_gpu_positions_lattice
)
