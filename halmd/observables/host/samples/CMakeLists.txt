halmd_add_library(halmd_observables_host_samples
  binned_phase_space.cpp
  phase_space.cpp
)
halmd_add_modules(
  libhalmd_observables_host_samples_binned_phase_space
  libhalmd_observables_host_samples_phase_space
)
