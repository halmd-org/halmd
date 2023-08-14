if(HALMD_WITH_GPU)
  set(DISABLE_GPU "--disable-gpu")
endif()

halmd_add_test(lua/mdsim/positions/excluded_volume
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/excluded_volume.lua
)
