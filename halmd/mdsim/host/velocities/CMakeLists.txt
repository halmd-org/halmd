halmd_add_library(halmd_mdsim_host_velocities
  boltzmann.cpp
  phase_space.cpp
)
halmd_add_modules(
  libhalmd_mdsim_host_velocities_boltzmann
  libhalmd_mdsim_host_velocities_phase_space
)
