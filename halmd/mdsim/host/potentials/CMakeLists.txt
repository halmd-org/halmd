halmd_add_library(
  halmd_mdsim_host_potentials_lennard_jones
  lennard_jones.cpp
)
halmd_add_modules(libhalmd_mdsim_host_potentials_lennard_jones)

halmd_add_potential(
  halmd_mdsim_host_potentials_lennard_jones_linear
  lennard_jones_linear
  lennard_jones_linear.cpp
)

halmd_add_potential(
  halmd_mdsim_host_potentials_modified_lennard_jones
  modified_lennard_jones
  modified_lennard_jones.cpp
)

halmd_add_potential(
  halmd_mdsim_host_potentials_morse
  morse
  morse.cpp
)

halmd_add_potential(
  halmd_mdsim_host_potentials_power_law
  power_law
  power_law.cpp
)

halmd_add_potential(
  halmd_mdsim_host_potentials_power_law_with_core
  power_law_with_core
  power_law_with_core.cpp
)
