add_executable(nvlock
  wrapper.c
  )
install(TARGETS nvlock
  RUNTIME DESTINATION bin
  )

add_library(libnvlock SHARED
  nvlock.cpp
  )
target_link_libraries(libnvlock
  ${CUDA_LIBRARY}
  )
set_target_properties(libnvlock PROPERTIES
  PREFIX ""
  )
install(TARGETS libnvlock
  LIBRARY DESTINATION lib
  )
