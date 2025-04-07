if(HALMD_WITH_GPU)
  set(DISABLE_GPU "--disable-gpu")
endif()

halmd_add_test(lua/recipes/create_mixture
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/create_mixture.lua --verbose
)

halmd_add_test(lua/recipes/create_slit_pore
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/create_slit_pore.lua --verbose
)
