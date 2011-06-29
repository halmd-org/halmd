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
  add_custom_target(halmd_doc_${OUTPUT_TYPE} ALL
    ${SPHINX_EXECUTABLE}
    -c .
    -b ${OUTPUT_TYPE}
    -d ${DOC_CACHE_DIRECTORY}
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${OUTPUT_TYPE}
    COMMENT "Building ${OUTPUT_TYPE} documentation"
    )
endforeach(OUTPUT_TYPE)

# serialise targets with a common cache directory
add_dependencies(halmd_doc_latex halmd_doc_html)

if(PDFLATEX_COMPILER)
  add_custom_target(halmd_doc_pdf ALL
    ${CMAKE_MAKE_PROGRAM}
    -C ${CMAKE_CURRENT_BINARY_DIR}/latex
    COMMENT "Building pdf documentation"
    )
  add_dependencies(halmd_doc_pdf halmd_doc_latex)
endif(PDFLATEX_COMPILER)

install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html
  DESTINATION share/doc/halmd
  PATTERN "cmake" EXCLUDE
  )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/latex/HALMD.pdf
  DESTINATION share/doc/halmd
  )

