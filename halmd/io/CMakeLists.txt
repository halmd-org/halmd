add_library(halmd_io
  logger.cpp
  )
target_link_libraries(halmd_io
  ${Boost_LOG_LIBRARY}
  ${Boost_THREAD_LIBRARY}
  )

add_subdirectory(profile)
add_subdirectory(trajectory)
add_subdirectory(statevars)
