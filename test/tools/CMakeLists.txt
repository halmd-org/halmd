add_library(test_tools_lua
  lua.cpp
)
target_link_libraries(test_tools_lua
  ${LUABIND_LIBRARIES}
  ${LUA_LIBRARIES}
  dl
)

# This test library contains an __attribute__((constructor)) function, which
# is not called explicitly by the executable it is linked against. To prevent
# the linker from optimizing away the function, we create a shared library.
add_library(test_tools_ctest SHARED
  ctest.cpp
)
