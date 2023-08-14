# Lua tests for module wavevector
halmd_add_test(unit/observables/utility/wavevector
  ${HALMD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/wavevector.lua --verbose --verbose
)

