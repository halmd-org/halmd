add_executable(test_thermodynamics
  thermodynamics.cpp
)
target_link_libraries(test_thermodynamics
  ${HALMD_COMMON_LIBRARIES}
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(thermodynamics test_thermodynamics)
