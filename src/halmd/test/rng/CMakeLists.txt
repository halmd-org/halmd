if(WITH_CUDA)
  add_executable(rand48_test rand48.cpp)
  target_link_libraries(rand48_test
    gpu_rand48
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${CUDA_LIBRARIES}
    pthread
    rt
    )
  foreach(count
    100
    1000
    10000
    100000
    1000000
    10000000
    100000000
  )
    add_test(rand48_test_${count} rand48_test -N${count})
  endforeach(count)
endif(WITH_CUDA)
