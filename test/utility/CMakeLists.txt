add_executable(test_demangle
  demangle.cpp
)
target_link_libraries(test_demangle
  ${HALMD_TEST_LIBRARIES}
)
add_test(demangle
  test_demangle --log_level=test_suite
)

add_executable(test_signal
  signal.cpp
)
target_link_libraries(test_signal
  ${HALMD_TEST_LIBRARIES}
)
add_test(signal
  test_signal --log_level=test_suite
)

add_executable(test_options_parser
  options_parser.cpp
)
target_link_libraries(test_options_parser
  test_tools_lua
  halmd_utility
  halmd_utility_program_options
  halmd_utility_lua_wrapper
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  ${HALMD_TEST_LIBRARIES}
)
add_test(options_parser
  test_options_parser --log_level=test_suite
)

add_executable(test_posix_signal
  posix_signal.cpp
)
target_link_libraries(test_posix_signal
  halmd_utility
  ${Boost_SYSTEM_LIBRARY}
  ${HALMD_TEST_LIBRARIES}
  rt
)
add_test(posix_signal
  test_posix_signal --log_level=message
)

add_subdirectory(lua_wrapper)
add_subdirectory(program_options)

if(WITH_CUDA)
  add_subdirectory(gpu)
endif(WITH_CUDA)
