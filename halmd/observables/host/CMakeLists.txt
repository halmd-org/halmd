halmd_add_library(halmd_observables_host
  density_mode.cpp
  current_density_mode.cpp
  phase_space.cpp
  thermodynamics.cpp
)
halmd_add_modules(
  libhalmd_observables_host_density_mode
  libhalmd_observables_host_current_density_mode
  libhalmd_observables_host_phase_space
  libhalmd_observables_host_thermodynamics
)

add_subdirectory(dynamics)
add_subdirectory(samples)
