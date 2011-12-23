##
# HALMD integration tests
#
add_subdirectory(continuation)
add_subdirectory(examples)

configure_file(libhalmd.lua.in libhalmd.lua @ONLY)
add_test(integration/libhalmd ${LUA_EXECUTABLE} libhalmd.lua)
