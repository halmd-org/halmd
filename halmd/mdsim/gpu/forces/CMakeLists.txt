add_library(halmd_mdsim_gpu_forces
  lennard_jones.cpp
  morse.cpp
  pair_trunc_potentials.cu
  zero.cpp
)
