add_executable(nvlock
  wrapper.c
  )
install(TARGETS nvlock
  RUNTIME DESTINATION bin
  )

add_library(libnvlock SHARED
  nvlock.c
  )
target_link_libraries(libnvlock
  dl
  )
set_target_properties(libnvlock PROPERTIES
  PREFIX ""
  )
install(TARGETS libnvlock
  LIBRARY DESTINATION lib
  )
