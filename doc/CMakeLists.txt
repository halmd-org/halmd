set(DOC_OUTPUT_TYPES
  html
  latex
  )
set(DOC_CACHE_DIRECTORY ".doctrees")

configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/_conf.py.in
  ${CMAKE_CURRENT_BINARY_DIR}/conf.py
  ESCAPE_QUOTES @ONLY
  )
set(DOC_OUTPUT_FILES
  ${DOC_OUTPUT_TYPES}
  ${DOC_CACHE_DIRECTORY}
  )
set_directory_properties(PROPERTIES
  ADDITIONAL_MAKE_CLEAN_FILES "${DOC_OUTPUT_FILES}"
  )
foreach(OUTPUT_TYPE ${DOC_OUTPUT_TYPES})
  add_custom_target(doc_${OUTPUT_TYPE} ALL
    ${SPHINX_EXECUTABLE}
    -c .
    -b ${OUTPUT_TYPE}
    -d ${DOC_CACHE_DIRECTORY}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${OUTPUT_TYPE}
    )
endforeach(OUTPUT_TYPE)
