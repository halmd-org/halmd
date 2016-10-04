halmd_add_potential(
  halmd_mdsim_host_potentials_pair_lennard_jones
  pair lennard_jones
  lennard_jones.cpp
)

halmd_add_potential(
  halmd_mdsim_host_potentials_pair_modified_lennard_jones
  pair modified_lennard_jones
  modified_lennard_jones.cpp
)

halmd_add_potential(
  halmd_mdsim_host_potentials_pair_morse
  pair morse
  morse.cpp
)

halmd_add_potential(
  halmd_mdsim_host_potentials_pair_power_law
  pair power_law
  power_law.cpp
)

halmd_add_potential(
  halmd_mdsim_host_potentials_pair_power_law_with_core
  pair power_law_with_core
  power_law_with_core.cpp
)
