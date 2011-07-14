add_subdirectory(integrators)
add_subdirectory(positions)
add_subdirectory(velocities)

add_executable(test_unit_mdsim_clock
  clock.cpp
)
target_link_libraries(test_unit_mdsim_clock
  libhalmd
  ${HALMD_COMMON_LIBRARIES}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/mdsim/clock
  test_unit_mdsim_clock --log_level=test_suite
)
