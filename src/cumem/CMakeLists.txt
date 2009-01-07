add_executable(cumem cumem.c)
target_link_libraries(cumem
  ${CUDA_LIBRARY}
  )
