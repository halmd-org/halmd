halmd_add_library(halmd_observables_host
  chemical_potential.cpp
  density_mode.cpp
  phase_space.cpp
  thermodynamics.cpp
)
halmd_add_modules(
  libhalmd_observables_host_chemical_potential
  libhalmd_observables_host_density_mode
  libhalmd_observables_host_phase_space
  libhalmd_observables_host_thermodynamics
)

add_subdirectory(dynamics)
add_subdirectory(samples)
