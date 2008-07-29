INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src/mdsim)

ADD_EXECUTABLE(sample_smooth_function smooth.cpp)
SET_TARGET_PROPERTIES(sample_smooth_function PROPERTIES COMPILE_FLAGS "-DDIM_3D -DUSE_SMOOTH_POTENTIAL")
TARGET_LINK_LIBRARIES(sample_smooth_function ljfluid3d_smooth_gpu ${Boost_PROGRAM_OPTIONS_LIBRARY} ${CUDA_LIBRARIES})

ADD_EXECUTABLE(parallel_scan scan.cpp)
TARGET_LINK_LIBRARIES(parallel_scan scan_gpu ${Boost_PROGRAM_OPTIONS_LIBRARY} ${CUDA_LIBRARIES})

ADD_EXECUTABLE(parallel_radix_sort radix.cpp)
TARGET_LINK_LIBRARIES(parallel_radix_sort radix_gpu scan_gpu rand48_gpu ljfluid3d_gpu ${Boost_PROGRAM_OPTIONS_LIBRARY} ${CUDA_LIBRARIES})

ADD_EXECUTABLE(rand48_test rand48.cpp)
TARGET_LINK_LIBRARIES(rand48_test rand48_gpu ljfluid3d_gpu scan_gpu ${Boost_PROGRAM_OPTIONS_LIBRARY} ${CUDA_LIBRARIES})
