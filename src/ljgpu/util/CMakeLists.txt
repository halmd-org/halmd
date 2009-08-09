add_library(util
  H5xx.cpp
  log.cpp
  signal.cpp
  )
target_link_libraries(util
  boost_log
  )
