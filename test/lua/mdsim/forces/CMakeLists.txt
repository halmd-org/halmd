if(HALMD_WITH_GPU)
  set(DISABLE_GPU "--disable-gpu")

  halmd_add_test(lua/mdsim/gpu/forces/tabulated_external
    ${HALMD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tabulated_external.lua --verbose
  )
endif()

halmd_add_test(lua/mdsim/host/forces/tabulated_external
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/tabulated_external.lua --verbose
)
