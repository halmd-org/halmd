halmd_add_library(halmd_mdsim_gpu_velocities
  boltzmann.cpp
  boltzmann_kernel.cu
)
halmd_add_modules(
  libhalmd_mdsim_gpu_velocities_boltzmann
)
