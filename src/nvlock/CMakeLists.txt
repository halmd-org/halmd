add_executable(nvlock nvlock.c)
target_link_libraries(nvlock
  ${CUDA_LIBRARY}
  )
install(TARGETS nvlock
  RUNTIME DESTINATION bin
  )
