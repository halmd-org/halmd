halmd_add_library(halmd_observables_gpu_samples
  binned_phase_space.cpp
  phase_space.cpp
)
halmd_add_modules(
  libhalmd_observables_gpu_samples_binned_phase_space
  libhalmd_observables_gpu_samples_phase_space
)
