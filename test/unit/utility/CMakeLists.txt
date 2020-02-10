add_executable(test_unit_utility_demangle
  demangle.cpp
)
target_link_libraries(test_unit_utility_demangle
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/utility/demangle
  test_unit_utility_demangle --log_level=test_suite
)

add_executable(test_unit_utility_cache
  cache.cpp
)
target_link_libraries(test_unit_utility_cache
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/utility/cache
  test_unit_utility_cache --log_level=test_suite
)

add_executable(test_unit_utility_multi_index
  multi_index.cpp
)
target_link_libraries(test_unit_utility_multi_index
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/utility/multi_index
  test_unit_utility_multi_index --log_level=test_suite
)

add_executable(test_unit_utility_profiler
  profiler.cpp
)
target_link_libraries(test_unit_utility_profiler
  halmd_utility
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/utility/profiler
  test_unit_utility_profiler --log_level=test_suite
)

add_executable(test_unit_utility_signal
  signal.cpp
)
target_link_libraries(test_unit_utility_signal
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/utility/signal
  test_unit_utility_signal --log_level=test_suite
)
# If the signal unit tests fail, they might fail with *** glibc detected ***
# memory corruption errors, which would freeze the test program infinitely.
set_property(TEST unit/utility/signal
  PROPERTY TIMEOUT 5
)

if (NOT ${Boost_VERSION} LESS 1.59.0)
  add_executable(test_unit_utility_posix_signal
    posix_signal.cpp
  )
  target_link_libraries(test_unit_utility_posix_signal
    halmd_utility
    ${HALMD_TEST_LIBRARIES}
  )
  add_test(unit/utility/posix_signal
    test_unit_utility_posix_signal --log_level=test_suite
  )
endif()

if (NOT ${Boost_VERSION} LESS 1.59.0)
  add_executable(test_unit_utility_raw_array
    raw_array.cpp
  )
  target_link_libraries(test_unit_utility_raw_array
    ${HALMD_TEST_LIBRARIES}
  )
  add_test(unit/utility/raw_array
    test_unit_utility_raw_array --log_level=test_suite
  )
endif()

add_subdirectory(lua)
if(HALMD_WITH_GPU)
  add_subdirectory(gpu)
endif()
