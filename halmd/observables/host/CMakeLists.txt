halmd_add_library(halmd_observables_host
  binned_phase_space.cpp
  density_mode.cpp
  phase_space.cpp
  profiles.cpp
  thermodynamics.cpp
)
halmd_add_modules(
  libhalmd_observables_host_binned_phase_space
  libhalmd_observables_host_density_mode
  libhalmd_observables_host_phase_space
  libhalmd_observables_host_profiles
  libhalmd_observables_host_thermodynamics
)

add_subdirectory(dynamics)
add_subdirectory(samples)
