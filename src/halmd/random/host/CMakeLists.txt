add_library(halmd_random_host
  random.cpp
)
target_link_libraries(halmd_random_host
  ${GSL_LIBRARIES}
)
