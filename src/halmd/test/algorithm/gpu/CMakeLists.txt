add_executable(test_scan
  scan_kernel.cu
  scan.cpp
)
target_link_libraries(test_scan
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${CUDA_LIBRARIES}
)
foreach(count 10000000 1000000 100000 10000 1024 1000 512 100 10 2)
  add_test(test_scan_${count}
    test_scan -N${count}
  )
endforeach(count)

add_executable(test_radix_sort
  radix_sort.cpp
)
target_link_libraries(test_radix_sort
  halmd_algorithm_gpu
  halmd_random_gpu
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${CUDA_LIBRARIES}
)
