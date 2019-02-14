configure_file(
  two_particles.cmake.in
  test_integration_modularity_two_particles.cmake
  @ONLY
)

add_test(integration/modularity/two_particles/host/2d
  ${CMAKE_COMMAND} -DDIMENSION=2 -DUSE_HOST=1 -P test_integration_modularity_two_particles.cmake
)
add_test(integration/modularity/two_particles/host/3d
  ${CMAKE_COMMAND} -DDIMENSION=3 -DUSE_HOST=1 -P test_integration_modularity_two_particles.cmake
)
if(HALMD_WITH_GPU)
  halmd_add_gpu_test(integration/modularity/two_particles/gpu/2d
    ${CMAKE_COMMAND} -DDIMENSION=2 -P test_integration_modularity_two_particles.cmake
  )
  halmd_add_gpu_test(integration/modularity/two_particles/gpu/3d
    ${CMAKE_COMMAND} -DDIMENSION=3 -P test_integration_modularity_two_particles.cmake
  )
endif()
