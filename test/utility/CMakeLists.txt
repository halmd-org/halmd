add_executable(test_demangle
  demangle.cpp
)
target_link_libraries(test_demangle
  ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY}
)
add_test(demangle
  test_demangle --log_level=test_suite
)
add_subdirectory(lua_wrapper)
add_subdirectory(program_options)
