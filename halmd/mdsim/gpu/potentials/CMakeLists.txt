halmd_add_library(
  halmd_mdsim_gpu_potentials_lennard_jones
  lennard_jones.cpp
  lennard_jones_kernel.cu
  lennard_jones_simple.cpp
  lennard_jones_simple_kernel.cu
)
halmd_add_modules(libhalmd_mdsim_gpu_potentials_lennard_jones)

halmd_add_potential(
  halmd_mdsim_gpu_potentials_lennard_jones_linear
  lennard_jones_linear
  lennard_jones_linear.cpp
  lennard_jones_linear_kernel.cu
)

halmd_add_potential(
  halmd_mdsim_gpu_potentials_modified_lennard_jones
  modified_lennard_jones
  modified_lennard_jones.cpp
  modified_lennard_jones_kernel.cu
)

halmd_add_potential(
  halmd_mdsim_gpu_potentials_morse
  morse
  morse.cpp
  morse_kernel.cu
)

halmd_add_potential(
  halmd_mdsim_gpu_potentials_power_law
  power_law
  power_law.cpp
  power_law_kernel.cu
)

halmd_add_potential(
  halmd_mdsim_gpu_potentials_power_law_with_core
  power_law_with_core
  power_law_with_core.cpp
  power_law_with_core_kernel.cu
)
