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

if(NOT GENERATE_DOC)
  message(WARNING "Missing required documentation tools")
  return()
endif()

configure_file(conf.py.in conf.py)

add_custom_target(halmd_doc_html ALL
  ${SPHINX_EXECUTABLE}
    -q -b html
    -c "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/html"
  COMMENT "Build HALMD HTML documentation"
)

configure_file(doxygen.cfg.in doxygen.cfg)

add_custom_target(halmd_doc_reference ALL
  COMMAND "${DOXYGEN}" doxygen.cfg
  COMMENT "Build HALMD C++ API reference"
  DEPENDS halmd_doc_html
)

install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html"
  DESTINATION "share/doc/halmd"
)

add_custom_target(halmd_doc_man ALL
  ${SPHINX_EXECUTABLE}
    -q -b man
    -c "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/man"
  COMMENT "Build HALMD man page"
  DEPENDS halmd_doc_html
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/man/halmd.1"
  DESTINATION "share/man/man1"
)

add_custom_target(halmd_doc_pdf ALL
  ${SPHINX_EXECUTABLE}
    -q -b latex
    -c "${CMAKE_CURRENT_BINARY_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/pdf"
  COMMENT "Build HALMD PDF documentation"
  DEPENDS halmd_doc_man
)
add_custom_command(TARGET halmd_doc_pdf POST_BUILD
  COMMAND ${CMAKE_MAKE_PROGRAM} LATEXOPTS=-interaction=batchmode
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/pdf"
)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/pdf/halmd.pdf"
  DESTINATION "share/doc/halmd"
)

set_directory_properties(PROPERTIES
  ADDITIONAL_MAKE_CLEAN_FILES ".buildinfo;.doctrees;html;man;pdf"
)
