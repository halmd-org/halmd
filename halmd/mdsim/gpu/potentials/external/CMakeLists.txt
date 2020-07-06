halmd_add_potential(
  halmd_mdsim_gpu_potentials_external_harmonic
  external harmonic
  harmonic.cpp
  harmonic_kernel.cu
)

halmd_add_potential(
  halmd_mdsim_gpu_potentials_external_planar_wall
  external planar_wall
  planar_wall.cpp
  planar_wall_kernel.cu
)
