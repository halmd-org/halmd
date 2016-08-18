halmd_add_library(halmd_observables_host_samples
  phase_space.cpp
  sample.cpp
)
halmd_add_modules(
  libhalmd_observables_host_samples_phase_space
  libhalmd_observables_host_samples_sample
)
