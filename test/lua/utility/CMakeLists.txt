if(HALMD_WITH_GPU)
  set(DISABLE_GPU "--disable-gpu")
endif()

halmd_add_test(lua/utility/string
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/string.lua --verbose
)

halmd_add_test(lua/utility/table
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/table.lua --verbose
)

