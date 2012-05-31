add_test(test/lua/mdsim/positions/excluded_volume/host
  ${HALMD_EXECUTABLE} --disable-gpu ${CMAKE_CURRENT_SOURCE_DIR}/excluded_volume.lua -vv
)
add_test(test/lua/mdsim/positions/excluded_volume/gpu
  ${HALMD_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/excluded_volume.lua -vv
)
