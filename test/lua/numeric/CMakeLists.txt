if(HALMD_WITH_GPU)
  set(DISABLE_GPU "--disable-gpu")
endif()

add_test(lua/numeric/multi_index
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/multi_index.lua
)

