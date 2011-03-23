add_executable(test_accumulator
  accumulator.cpp
)
target_link_libraries(test_accumulator
  ${HALMD_TEST_LIBRARIES}
)
add_test(numeric_accumulator
  test_accumulator --log_level=test_suite
)

add_executable(test_pow
  pow.cpp
)
target_link_libraries(test_pow
  ${HALMD_TEST_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  rt
  )
add_test(numeric_pow
  test_pow --log_level=test_suite
)

add_subdirectory(mp)
