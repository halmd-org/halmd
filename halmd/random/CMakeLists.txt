halmd_add_module(halmd_random_random
  random.cpp
)

add_subdirectory(host)
if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)

