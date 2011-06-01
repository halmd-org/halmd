# Generate documentation in HTML and PDF format using Sphinx.

set(GENERATE_DOC TRUE)

# We use the Sphinx documentation generator to render HTML and manual
# pages from the user and reference documentation in ReST format.
find_package(Sphinx QUIET)
if(NOT SPHINX_FOUND)
  message(WARNING "Unable to find Sphinx documentation generator")
  set(GENERATE_DOC FALSE)
endif(NOT SPHINX_FOUND)

# We use Doxygen to extract inline reference documentation from the
# source code in XML format. This command also finds the dot program
# of the Graphviz package to render class inheritance diagrams.
find_package(Doxygen QUIET)
if(NOT DOXYGEN_FOUND)
  message(WARNING "Unable to find Doxygen documentation generator")
  set(GENERATE_DOC FALSE)
endif(NOT DOXYGEN_FOUND)

# The XML output is transformed to ReST using XSLT.
find_program(XSLTPROC_EXECUTABLE xsltproc)
if(NOT XSLTPROC_EXECUTABLE)
  message(WARNING "Unable to find xsltproc command")
  set(GENERATE_DOC FALSE)
endif(NOT XSLTPROC_EXECUTABLE)
mark_as_advanced(XSLTPROC_EXECUTABLE)

# Sphinx uses latex and dvipng to render LaTeX formulas.
find_package(LATEX QUIET)
if(NOT LATEX_COMPILER)
  message(WARNING "Unable to find LaTeX compiler")
  set(GENERATE_DOC FALSE)
endif(NOT LATEX_COMPILER)
find_program(DVIPNG_CONVERTER dvipng)
if(NOT DVIPNG_CONVERTER)
  message(WARNING "Unable to find dvipng converter")
  set(GENERATE_DOC FALSE)
endif(NOT DVIPNG_CONVERTER)
mark_as_advanced(DVIPNG_CONVERTER)

if(GENERATE_DOC)
  # documentation tools
  set(SOURCE_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/_build")
  # configured documentation tools and intermediate build results
  set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")
  # static ReST documentation sources
  set(SOURCES_DIR "${CMAKE_CURRENT_BINARY_DIR}/_sources")
  # generated ReST documentation sources
  set(REF_SOURCES_DIR "${SOURCES_DIR}/reference")
  # master document with modules index
  set(REF_MASTER_DOC "modules")

  # substitute variables in configuration and scripts
  foreach(file
      conf.py
      reference.dox
      reference.xsl
      sources.cmake
  )
    configure_file(
      "${SOURCE_BUILD_DIR}/${file}.in"
      "${BINARY_BUILD_DIR}/${file}"
      @ONLY
    )
  endforeach(file)

  add_custom_command(
    OUTPUT "${BINARY_BUILD_DIR}/reference.xml"
    COMMAND "${DOXYGEN}" reference.dox
    COMMAND "${XSLTPROC_EXECUTABLE}" -o reference.xml xml/combine.xslt xml/index.xml
    WORKING_DIRECTORY "${BINARY_BUILD_DIR}"
    DEPENDS "${BINARY_BUILD_DIR}/reference.dox"
    COMMENT "Extracting inline documentation"
  )
  set(CLEAN_FILES
    "${BINARY_BUILD_DIR}/xml"
    "${BINARY_BUILD_DIR}/html"
  )
  add_custom_command(
    OUTPUT "${REF_SOURCES_DIR}/${REF_MASTER_DOC}.rst"
    COMMAND "${XSLTPROC_EXECUTABLE}" -o "${REF_SOURCES_DIR}/" reference.xsl reference.xml
    WORKING_DIRECTORY "${BINARY_BUILD_DIR}"
    DEPENDS "${BINARY_BUILD_DIR}/reference.xsl" "${BINARY_BUILD_DIR}/reference.xml"
    COMMENT "Building reference documentation"
  )
  add_custom_target(halmd_doc_reference ALL
    DEPENDS "${REF_SOURCES_DIR}/${REF_MASTER_DOC}.rst"
  )

  # Sphinx requires all sources in the same directory tree. As we wish
  # to include generated reference documention from the build tree, we
  # copy static ReST documents to the build tree before calling Sphinx.
  add_custom_target(halmd_doc_sources ALL
    "${CMAKE_COMMAND}" -P "${BINARY_BUILD_DIR}/sources.cmake"
  )
  list(APPEND CLEAN_FILES
    "${SOURCES_DIR}"
  )
  add_dependencies(halmd_doc_sources
    halmd_doc_reference
  )
  # note the trailing slash to exclude directory name
  install(DIRECTORY "${SOURCES_DIR}/"
    DESTINATION "share/doc/halmd"
  )

  # Sphinx cache with pickled ReST documents
  set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")
  # HTML output directory
  set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")

  # This target builds HTML documentation using Sphinx.
  add_custom_target(halmd_doc_html ALL
    ${SPHINX_EXECUTABLE}
      -q -b html
      -c "${BINARY_BUILD_DIR}"
      -d "${SPHINX_CACHE_DIR}"
      "${SOURCES_DIR}"
      "${SPHINX_HTML_DIR}"
    COMMENT "Building HTML documentation with Sphinx"
  )
  list(APPEND CLEAN_FILES
    "${SPHINX_CACHE_DIR}"
    "${SPHINX_HTML_DIR}"
  )
  add_dependencies(halmd_doc_html
    halmd_doc_sources
  )
  install(DIRECTORY "${SPHINX_HTML_DIR}"
    DESTINATION "share/doc/halmd"
  )

  # HTML output directory
  set(SPHINX_MAN_DIR "${CMAKE_CURRENT_BINARY_DIR}/man")
  # This target builds a manual page using Sphinx.

  add_custom_target(halmd_doc_man ALL
    ${SPHINX_EXECUTABLE}
      -q -b man
      -c "${BINARY_BUILD_DIR}"
      -d "${SPHINX_CACHE_DIR}"
      "${SOURCES_DIR}"
      "${SPHINX_MAN_DIR}"
    COMMENT "Building manual page with Sphinx"
  )
  list(APPEND CLEAN_FILES
    "${SPHINX_MAN_DIR}"
  )
  add_dependencies(halmd_doc_man
    halmd_doc_sources
  )
  # serialize Sphinx targets to avoid cache conflicts in parallel builds
  add_dependencies(halmd_doc_man
    halmd_doc_html
  )
  install(FILES "${SPHINX_MAN_DIR}/halmd.1"
    DESTINATION "share/man/man1"
  )

  # This target builds PDF documentation using Sphinx and LaTeX.
  if(PDFLATEX_COMPILER)
    # PDF output directory
    set(SPHINX_PDF_DIR "${CMAKE_CURRENT_BINARY_DIR}/pdf")

    add_custom_target(halmd_doc_pdf ALL
      ${SPHINX_EXECUTABLE}
	-q -b latex
	-c "${BINARY_BUILD_DIR}"
	-d "${SPHINX_CACHE_DIR}"
	"${SOURCES_DIR}"
	"${SPHINX_PDF_DIR}"
      COMMENT "Building PDF documentation with Sphinx"
    )
    add_custom_command(TARGET halmd_doc_pdf POST_BUILD
      COMMAND ${CMAKE_MAKE_PROGRAM} LATEXOPTS=-interaction=batchmode
      WORKING_DIRECTORY "${SPHINX_PDF_DIR}"
    )
    list(APPEND CLEAN_FILES
      "${SPHINX_PDF_DIR}"
    )
    add_dependencies(halmd_doc_pdf
      halmd_doc_sources
    )
    # serialize Sphinx targets to avoid cache conflicts in parallel builds
    add_dependencies(halmd_doc_pdf
      halmd_doc_man
    )
    install(FILES "${SPHINX_PDF_DIR}/halmd.pdf"
      DESTINATION "share/doc/halmd"
    )
  endif(PDFLATEX_COMPILER)

  # Add output directories to clean target.
  set_directory_properties(PROPERTIES
    ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}"
  )

else(GENERATE_DOC)
  message(WARNING "Missing required documentation tools")
endif(GENERATE_DOC)
