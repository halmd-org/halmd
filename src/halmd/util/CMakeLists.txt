add_library(util
  H5xx.cpp
  log.cpp
  signal.cpp
  )
target_link_libraries(util
  ${Boost_LOG_LIBRARY}
  )
