halmd_add_library(halmd_observables
  runtime_estimate.cpp
  sampler.cpp
  ssf.cpp
  thermodynamics.cpp
)
halmd_add_modules(
  libhalmd_observables_runtime_estimate
  libhalmd_observables_sampler
  libhalmd_observables_ssf
  libhalmd_observables_thermodynamics
)

add_subdirectory(samples) # must appear first in order to register blocking_scheme_base before blocking_scheme
add_subdirectory(dynamics)
add_subdirectory(host)
if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)
add_subdirectory(utility)
