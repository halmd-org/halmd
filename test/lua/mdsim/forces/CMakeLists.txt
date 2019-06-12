add_test(lua/mdsim/host/forces/tabulated_external
  ${HALMD_EXECUTABLE} --disable-gpu ${CMAKE_CURRENT_SOURCE_DIR}/tabulated_external.lua
)
if(HALMD_WITH_GPU)
add_test(lua/mdsim/gpu/forces/tabulated_external
  ${HALMD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tabulated_external.lua
)
endif()
