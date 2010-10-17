add_library(test_tools_lua
  lua.cpp
)
target_link_libraries(test_tools_lua
  ${LUABIND_LIBRARIES}
  ${LUA_LIBRARIES}
  dl
)
