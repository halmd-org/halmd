##
# HALMD Lua unit tests
#
add_subdirectory(mdsim)
add_subdirectory(numeric)
add_subdirectory(recipes)

# Lua test for simulation script
configure_file(
  script.cmake.in
  test_lua_script.cmake
  @ONLY
)

add_test(lua/script
  ${CMAKE_COMMAND} -P test_lua_script.cmake
)
