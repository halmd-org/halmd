halmd_add_potential(
  halmd_mdsim_gpu_potentials_pair_lennard_jones
  pair lennard_jones
  lennard_jones.cpp
  lennard_jones_kernel.cu
  lennard_jones_simple.cpp
  lennard_jones_simple_kernel.cu
)
if(HALMD_WITH_pair_lennard_jones)
  # the "simple" version needs to be loaded separately
  halmd_add_modules("libhalmd_mdsim_gpu_potentials_pair_lennard_jones_simple")
endif()

halmd_add_potential(
  halmd_mdsim_gpu_potentials_pair_lennard_jones_linear
  pair lennard_jones_linear
  lennard_jones_linear.cpp
  lennard_jones_linear_kernel.cu
)

halmd_add_potential(
  halmd_mdsim_gpu_potentials_pair_modified_lennard_jones
  pair modified_lennard_jones
  modified_lennard_jones.cpp
  modified_lennard_jones_kernel.cu
)

halmd_add_potential(
  halmd_mdsim_gpu_potentials_pair_morse
  pair morse
  morse.cpp
  morse_kernel.cu
)

halmd_add_potential(
  halmd_mdsim_gpu_potentials_pair_power_law
  pair power_law
  power_law.cpp
  power_law_kernel.cu
)

halmd_add_potential(
  halmd_mdsim_gpu_potentials_pair_power_law_with_core
  pair power_law_with_core
  power_law_with_core.cpp
  power_law_with_core_kernel.cu
)
