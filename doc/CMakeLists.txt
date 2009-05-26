if(DOXYGEN_FOUND)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/doxyfile
    ESCAPE_QUOTES @ONLY
    )
  add_custom_target(doc ALL
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/doxyfile
    )
  set(DOXYGEN_OUTPUT_FILES
    html
    latex
    )
  set_directory_properties(PROPERTIES
    ADDITIONAL_MAKE_CLEAN_FILES "${DOXYGEN_OUTPUT_FILES}"
    )
endif(DOXYGEN_FOUND)
