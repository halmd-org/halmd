if(HALMD_WITH_GPU)
  set(DISABLE_GPU "--disable-gpu")
endif()

add_test(lua/recipes/create_mixture
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/create_mixture.lua
)
