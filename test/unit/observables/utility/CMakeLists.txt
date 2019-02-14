# Lua tests for module wavevector
add_test(unit/observables/utility/wavevector
  ${HALMD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wavevector.lua --verbose --verbose
)

