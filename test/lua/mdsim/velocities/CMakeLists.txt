# FIXME add different number of particles

#
# module rescale
#
if(HALMD_WITH_GPU)
  halmd_add_test(lua/velocities/rescale/gpu/2d
    ${HALMD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/rescale.lua --verbose --dimension 2
  )
  halmd_add_test(lua/velocities/rescale/gpu/3d
    ${HALMD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/rescale.lua --verbose --dimension 3
  )
  set(DISABLE_GPU "--disable-gpu")
endif()

halmd_add_test(lua/velocities/rescale/host/2d
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/rescale.lua --verbose --dimension 2
)
halmd_add_test(lua/velocities/rescale/host/3d
  ${HALMD_EXECUTABLE} ${DISABLE_GPU} ${CMAKE_CURRENT_SOURCE_DIR}/rescale.lua --verbose --dimension 3
)
