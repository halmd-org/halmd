halmd_add_library(halmd_observables_host
  current_density_mode.cpp
  density_mode.cpp
  kinetic_energy_density_mode.cpp
  phase_space.cpp
  thermodynamics.cpp
)
halmd_add_modules(
  libhalmd_observables_host_current_density_mode
  libhalmd_observables_host_density_mode
  libhalmd_observables_host_kinetic_energy_density_mode
  libhalmd_observables_host_phase_space
  libhalmd_observables_host_thermodynamics
)

add_subdirectory(dynamics)
add_subdirectory(samples)
