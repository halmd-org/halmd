halmd_add_module(halmd_mdsim_gpu_velocities_boltzmann
  boltzmann.cpp
  boltzmann_kernel.cu
)
halmd_add_module(halmd_mdsim_gpu_velocities_phase_space
  phase_space.cpp
)
