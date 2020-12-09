if(HALMD_WITH_GPU)
  add_test(lua/particle_groups/region/gpu
    ${HALMD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/region.lua
  )
  set(DISABLE_GPU "--disable-gpu")
endif()

add_test(lua/particle_groups/region/host
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/region.lua
)

# FIXME add basic tests for all, id_range, region_species
