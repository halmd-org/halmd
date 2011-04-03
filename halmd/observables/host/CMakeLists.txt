halmd_add_library(halmd_observables_host
  mean_quartic_displacement.cpp
  mean_square_displacement.cpp
  density_mode.cpp
  phase_space.cpp
  thermodynamics.cpp
  velocity_autocorrelation.cpp
)
halmd_add_modules(
  libhalmd_observables_host_mean_quartic_displacement
  libhalmd_observables_host_mean_square_displacement
  libhalmd_observables_host_density_mode
  libhalmd_observables_host_phase_space
  libhalmd_observables_host_thermodynamics
  libhalmd_observables_host_velocity_autocorrelation
)

add_subdirectory(samples)
