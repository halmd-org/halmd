add_library(halmd_observables_host
  density_modes.cpp
  ssf.cpp
  thermodynamics.cpp
  trajectory.cpp
)
add_subdirectory(samples)
