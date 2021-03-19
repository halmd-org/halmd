# FIXME add basic tests for all, id_range

#
# module region
#
if(HALMD_WITH_GPU)
  add_test(lua/particle_groups/region/gpu/2d
    ${HALMD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/region.lua --verbose --dimension 2
  )
  add_test(lua/particle_groups/region/gpu/3d
    ${HALMD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/region.lua --verbose --dimension 3
  )
  set(DISABLE_GPU "--disable-gpu")
endif()

add_test(lua/particle_groups/region/host/2d
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/region.lua --verbose --dimension 2
)
add_test(lua/particle_groups/region/host/3d
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/region.lua --verbose --dimension 3
)

#
# module region species
#
if(HALMD_WITH_GPU)
  add_test(lua/particle_groups/region_species/gpu/2d
    ${HALMD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/region_species.lua --verbose --dimension 2
  )
  add_test(lua/particle_groups/region_species/gpu/3d
    ${HALMD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/region_species.lua --verbose --dimension 3
  )
  set(DISABLE_GPU "--disable-gpu")
endif()

add_test(lua/particle_groups/region_species/host/2d
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/region_species.lua --verbose --dimension 2
)
add_test(lua/particle_groups/region_species/host/3d
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/region_species.lua --verbose --dimension 3
)

