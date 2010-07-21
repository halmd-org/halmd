add_executable(test_accumulator
  accumulator.cpp
)
target_link_libraries(test_accumulator
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(accumulator
  test_accumulator
)

add_subdirectory(host)
add_subdirectory(mp)
