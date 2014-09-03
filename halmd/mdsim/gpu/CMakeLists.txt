halmd_add_modules(
  libhalmd_mdsim_gpu_binning
  libhalmd_mdsim_gpu_max_displacement
  libhalmd_mdsim_gpu_neighbour
  libhalmd_mdsim_gpu_particle
  libhalmd_mdsim_gpu_particle_group
)

add_subdirectory(integrators)
add_subdirectory(neighbours)
add_subdirectory(particle_groups)
add_subdirectory(positions)
add_subdirectory(potentials)
add_subdirectory(sorts)
add_subdirectory(velocities)

halmd_add_library(halmd_mdsim_gpu
  binning.cpp
  binning_kernel.cu
  max_displacement.cpp
  max_displacement_kernel.cu
  neighbour.cpp
  particle.cpp
  particle_group.cpp
  particle_kernel.cu
  velocity_kernel.cu
  particle_group_kernel.cu
)
