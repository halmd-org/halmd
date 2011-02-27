add_library(halmd_mdsim_gpu_velocities
  boltzmann.cpp
  boltzmann_kernel.cu
  trajectory.cpp
)
