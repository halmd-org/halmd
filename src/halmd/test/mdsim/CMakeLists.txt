add_executable(test_thermodynamics
  thermodynamics.cpp
)
target_link_libraries(test_thermodynamics
  ${HALMD_COMMON_LIBRARIES}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(thermodynamics_host_2d
  test_thermodynamics --run_test=host/2d/* --log_level=test_suite
)
add_test(thermodynamics_host_3d
  test_thermodynamics --run_test=host/3d/* --log_level=test_suite
)
add_test(thermodynamics_gpu_2d
  test_thermodynamics --run_test=gpu/2d/* --log_level=test_suite
)
add_test(thermodynamics_gpu_3d
  test_thermodynamics --run_test=gpu/3d/* --log_level=test_suite
)
