# clean build tree subdirectory
file(GLOB_RECURSE halmd_lua_files "${CMAKE_CURRENT_BINARY_DIR}/halmd/*.lua")
foreach(file ${halmd_lua_files})
  file(REMOVE ${file})
endforeach(file)

# copy Lua source files from source to build tree
file(GLOB_RECURSE halmd_lua_sources RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "halmd/*.lua.in")
foreach(file ${halmd_lua_sources})
  string(REGEX REPLACE "\\.in$" "" out_file ${file})
  configure_file(${file} ${out_file} @ONLY)
endforeach(file)

# install Lua files from build tree
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/halmd"
  DESTINATION share
  FILES_MATCHING PATTERN "*.lua"
)
