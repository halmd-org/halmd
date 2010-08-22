set(_HALMD_GENERATE_DOC TRUE)

# We use the Sphinx documentation generator to render HTML and manual
# pages from the user and reference documentation in ReST format.
find_package(Sphinx QUIET)
if(NOT SPHINX_FOUND)
  message(STATUS "Unable to find Sphinx documentation generator")
  set(_HALMD_GENERATE_DOC FALSE)
endif(NOT SPHINX_FOUND)

# We use Doxygen to extract inline reference documentation from the
# source code in XML format. This command also finds the dot program
# of the Graphviz package to render class inheritance diagrams.
find_package(Doxygen QUIET)
if(NOT DOXYGEN_FOUND)
  message(STATUS "Unable to find Doxygen documentation generator")
  set(_HALMD_GENERATE_DOC FALSE)
endif(NOT DOXYGEN_FOUND)

# The XML output is transformed to ReST using XSLT.
find_program(XSLTPROC_EXECUTABLE xsltproc)
if(NOT XSLTPROC_EXECUTABLE)
  message(STATUS "Unable to find xsltproc command")
  set(_HALMD_GENERATE_DOC FALSE)
endif(NOT XSLTPROC_EXECUTABLE)
mark_as_advanced(XSLTPROC_EXECUTABLE)

# Sphinx uses latex and dvipng to render LaTeX formulas.
find_package(LATEX QUIET)
if(NOT LATEX_COMPILER)
  message(STATUS "Unable to find LaTeX compiler")
  set(_HALMD_GENERATE_DOC FALSE)
endif(NOT LATEX_COMPILER)
find_program(DVIPNG_CONVERTER dvipng)
if(NOT DVIPNG_CONVERTER)
  message(STATUS "Unable to find dvipng converter")
  set(_HALMD_GENERATE_DOC FALSE)
endif(NOT DVIPNG_CONVERTER)
mark_as_advanced(DVIPNG_CONVERTER)

if(NOT _HALMD_GENERATE_DOC)
  message(STATUS "Unable to generate ${PROJECT_NAME} documentation")
endif(NOT _HALMD_GENERATE_DOC)

# set global variable to indicate whether all tools were found
set(HALMD_GENERATE_DOC ${_HALMD_GENERATE_DOC}
  CACHE INTERNAL "Generate HALMD documentation" FORCE
)

if(HALMD_GENERATE_DOC)
  configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    ESCAPE_QUOTES @ONLY
    )
  add_custom_target(halmd_doc_xml ALL
    ${DOXYGEN_EXECUTABLE} Doxyfile
    COMMENT "Building doxygen xml output"
    )
  set(DOC_OUTPUT_FILES
    "xml"
    )

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
    ${DOC_OUTPUT_FILES}
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
    add_dependencies(halmd_doc_${OUTPUT_TYPE} halmd_doc_xml halmd_doc_impl)
  endforeach(OUTPUT_TYPE)

  add_custom_target(halmd_doc_impl
    COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/impl/get_modules.py
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/impl
    COMMENT "Extracting module structure and documentation"
    )

  # serialise targets with a common cache directory
  add_dependencies(halmd_doc_latex halmd_doc_html)

  if(PDFLATEX_COMPILER)
    add_custom_target(halmd_doc_pdf ALL
      ${CMAKE_MAKE_PROGRAM}
      -C ${CMAKE_CURRENT_BINARY_DIR}/latex
      LATEXOPTS=-interaction=batchmode
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
endif(HALMD_GENERATE_DOC)
