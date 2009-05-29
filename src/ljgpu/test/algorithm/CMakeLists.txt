if(WITH_CUDA)
  add_executable(parallel_scan scan.cpp)
  target_link_libraries(parallel_scan
    gpu_prefix_sum
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${CUDA_LIBRARIES}
    )
  foreach(count 10000000 1000000 100000 10000 1024 1000 512 100 10 2)
    add_test(parallel_scan_${count} parallel_scan -N${count})
  endforeach(count)

  add_executable(parallel_radix_sort radix.cpp)
  target_link_libraries(parallel_radix_sort
    gpu_radix_sort
    gpu_rand48
    ${Boost_PROGRAM_OPTIONS_LIBRARY}
    ${CUDA_LIBRARIES}
    )
endif(WITH_CUDA)
