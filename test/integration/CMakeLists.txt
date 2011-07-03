##
# HALMD integration tests
#
add_subdirectory(continuation)

configure_file(libhalmd.lua.in libhalmd.lua @ONLY)
add_test(integration/libhalmd lua libhalmd.lua)
