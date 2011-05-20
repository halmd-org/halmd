halmd_add_library(halmd_observables_dynamics
  blocking_scheme.cpp
  correlation.cpp
)
halmd_add_modules(
  libhalmd_observables_dynamics_blocking_scheme
  libhalmd_observables_dynamics_correlation
)
