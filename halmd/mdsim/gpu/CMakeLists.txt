add_library(halmd_mdsim_gpu
  force.cpp
  neighbour.cpp
  neighbour_kernel.cu
  particle.cpp
  particle_kernel.cu
  velocity.cpp
  velocity_kernel.cu
)

add_subdirectory(forces)
add_subdirectory(integrators)
add_subdirectory(positions)
add_subdirectory(sorts)
add_subdirectory(velocities)
