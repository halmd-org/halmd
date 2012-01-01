##
# HALMD integration tests
#
add_subdirectory(continuation)
add_subdirectory(examples)

add_test(integration/libhalmd
  ${LUA_EXECUTABLE} -lctest ${CMAKE_CURRENT_SOURCE_DIR}/libhalmd.lua
)
set_property(TEST integration/libhalmd
  PROPERTY ENVIRONMENT "LUA_PATH=${HALMD_LUA_PATH}" "LUA_CPATH=${HALMD_LUA_CPATH}"
)
