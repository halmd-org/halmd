if("${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}" VERSION_LESS "1.60")
  set(BOOST_TEST_ARGS "--log_level=test_suite")
else()
  set(BOOST_TEST_ARGS "--log_level=test_suite --")
endif()

add_subdirectory(continuation)
add_subdirectory(modularity)
add_subdirectory(validation)
