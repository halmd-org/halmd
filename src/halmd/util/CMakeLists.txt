add_library(util
  H5xx.cpp
  signal.cpp
  )
target_link_libraries(util
  ${Boost_THREAD_LIBRARY}
  )
