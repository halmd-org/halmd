add_library(halmd_rng
  random.cpp
  )
target_link_libraries(halmd_rng
  halmd_rng_host
  )

if(WITH_CUDA)
  target_link_libraries(halmd_rng
    halmd_rng_gpu
  )
endif(WITH_CUDA)

add_subdirectory(host)

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)
