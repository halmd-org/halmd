configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/Jamroot.in
  ${CMAKE_CURRENT_BINARY_DIR}/Jamroot
  ESCAPE_QUOTES @ONLY
  )
add_custom_target(halmd_doc ALL
  bjam --v2
  COMMENT "Building documentation using Boost.Jam"
  )
set_directory_properties(PROPERTIES
  ADDITIONAL_MAKE_CLEAN_FILES "bin html"
  )
install(DIRECTORY html
  DESTINATION share/doc/halmd
  )
