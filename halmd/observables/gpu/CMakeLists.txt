halmd_add_library(halmd_observables_gpu
  density_mode.cpp
  density_mode_kernel.cu
  phase_space.cpp
  phase_space_kernel.cu
  thermodynamics.cpp
  thermodynamics_kernel.cu
)
halmd_add_modules(
  libhalmd_observables_gpu_density_mode
  libhalmd_observables_gpu_phase_space
  libhalmd_observables_gpu_thermodynamics
)

add_subdirectory(dynamics)
add_subdirectory(samples)
