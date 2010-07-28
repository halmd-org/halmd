add_library(halmd_io
  logger.cpp
  )
target_link_libraries(halmd_io
  ${Boost_THREAD_LIBRARY}
  ${Boost_LOG_LIBRARY}
  )

add_subdirectory(profile)
add_subdirectory(trajectory)
