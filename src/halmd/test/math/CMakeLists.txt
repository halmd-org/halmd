add_executable(test_pow pow.cpp)
target_link_libraries(test_pow
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
  )
add_test("math::pow" test_pow)

add_subdirectory(dsfun)
