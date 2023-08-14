halmd_add_potential(
  halmd_mdsim_host_potentials_pair_custom
  pair custom
  custom.cpp
)

halmd_add_potential(
  halmd_mdsim_host_potentials_pair_lennard_jones
  pair lennard_jones
  lennard_jones.cpp
)

halmd_add_potential(
  halmd_mdsim_host_potentials_pair_mie
  pair mie
  mie.cpp
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
