halmd_add_library(halmd_random
  random.cpp
)
halmd_add_modules(
  libhalmd_random_random
)

add_subdirectory(host)
if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)
