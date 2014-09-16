halmd_add_potential(
  halmd_mdsim_gpu_potentials_external_harmonic
  external harmonic
  harmonic.cpp
  harmonic_kernel.cu
)

halmd_add_potential(
  halmd_mdsim_gpu_potentials_external_slit
  external slit
  slit.cpp
  slit_kernel.cu
)
