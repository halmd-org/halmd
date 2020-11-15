if(HALMD_WITH_GPU)
  set(DISABLE_GPU "--disable-gpu")
endif()

add_test(lua/utility/string
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/string.lua --verbose
)

