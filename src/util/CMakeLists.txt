add_executable(nvlock nvlock.c)
target_link_libraries(nvlock
  ${CUDA_LIBRARY}
  )
add_executable(cumem cumem.c)
target_link_libraries(cumem
  ${CUDA_LIBRARY}
  )

install(TARGETS nvlock
  RUNTIME DESTINATION bin
  )
