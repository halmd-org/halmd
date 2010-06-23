add_executable(test_thermodynamics
  thermodynamics.cpp
)
target_link_libraries(test_thermodynamics
  ${HALMD_COMMON_LIBRARIES}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(thermodynamics_host
  test_thermodynamics --run_test=host/* --log_level=test_suite
)
add_test(thermodynamics_gpu
  test_thermodynamics --run_test=gpu/* --log_level=test_suite
)
