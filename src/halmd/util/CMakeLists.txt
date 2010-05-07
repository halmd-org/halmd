add_library(util
  signal.cpp
  )
target_link_libraries(util
  ${Boost_THREAD_LIBRARY}
  )
