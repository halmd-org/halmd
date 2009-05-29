if(WITH_CUDA)
  add_executable(rand48_test rand48.cpp)
  target_link_libraries(rand48_test
    gpu_rand48
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${CUDA_LIBRARIES}
    )
endif(WITH_CUDA)
