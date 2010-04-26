add_library(halmd_utility
  options.cpp
  )
target_link_libraries(halmd_utility
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)
