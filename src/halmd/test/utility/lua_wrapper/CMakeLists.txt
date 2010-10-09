add_library(test_lua_setup
  lua_setup.cpp
)
target_link_libraries(test_lua_setup
  ${LUABIND_LIBRARIES}
  ${LUA_LIBRARIES}
)
