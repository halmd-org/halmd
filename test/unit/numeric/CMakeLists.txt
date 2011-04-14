add_executable(test_unit_numeric_accumulator
  accumulator.cpp
)
target_link_libraries(test_unit_numeric_accumulator
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/numeric/accumulator
  test_unit_numeric_accumulator --log_level=test_suite
)

add_executable(test_unit_numeric_pow
  pow.cpp
)
target_link_libraries(test_unit_numeric_pow
  ${HALMD_TEST_LIBRARIES}
  ${Boost_SYSTEM_LIBRARY}
  rt
  )
add_test(unit/numeric/pow
  test_unit_numeric_pow --log_level=test_suite
)

add_subdirectory(mp)
