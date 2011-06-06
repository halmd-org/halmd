add_executable(test_unit_utility_demangle
  demangle.cpp
)
target_link_libraries(test_unit_utility_demangle
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/utility/demangle
  test_unit_utility_demangle --log_level=message
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

add_executable(test_unit_utility_options_parser
  options_parser.cpp
)
target_link_libraries(test_unit_utility_options_parser
  test_tools_lua
  halmd_utility
  halmd_utility_program_options
  halmd_utility_lua
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${HALMD_TEST_LIBRARIES}
)
add_test(unit/utility/options_parser
  test_unit_utility_options_parser --log_level=test_suite
)

add_executable(test_unit_utility_posix_signal
  posix_signal.cpp
)
target_link_libraries(test_unit_utility_posix_signal
  halmd_utility
  ${Boost_SYSTEM_LIBRARY}
  ${HALMD_TEST_LIBRARIES}
  rt
)
add_test(unit/utility/posix_signal
  test_unit_utility_posix_signal --log_level=message
)

add_subdirectory(lua)
add_subdirectory(program_options)

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)