halmd_add_library(halmd_observables_dynamics
  correlation.cpp
  tcf.cpp
)
halmd_add_modules(
  libhalmd_observables_dynamics_correlation
  libhalmd_observables_dynamics_tcf_base
)
