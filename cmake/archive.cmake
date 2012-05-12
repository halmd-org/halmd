##
## Generate HALMD archive
##
# This set of rules crafts a tarball of the halmd repository, including git
# submodules and generated documentation as plain-text ReST, HTML with images,
# PDF and manual page.
#
# The tarball is created with git archive of the halmd repository, and further
# tarballs of submodules, and the documentation, are appended using GNU tar.
# This procedure ensures that the commit id of the HALMD commit is preserved,
# and may later be queried with
#
# bzip2 -d < halmd-<version>.tar.bz2 | git get-tar-commit-id
#
# If any command fails, the entire chain of commands is repeated,
# to ensure that all source and binary files are up-to-date.
#
# Caveats: Before creating an archive, ensure that all submodules are
# at the correct commit (as defined in the halmd commit), and purge
# any previously built documentation with make -C doc clean.
#

set(HALMD_ARCHIVE_PREFIX ${PROGRAM_NAME}-${HALMD_GIT_COMMIT_TAG})
set(HALMD_ARCHIVE_OUTPUT ${HALMD_ARCHIVE_PREFIX}.tar.bz2)

set(HALMD_ARCHIVE_INTERMEDIATE_OUTPUT CMakeFiles/archive_${PROGRAM_NAME})

add_custom_command(OUTPUT ${HALMD_ARCHIVE_INTERMEDIATE_OUTPUT}
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  COMMAND git archive
    --format=tar
    --output=${CMAKE_CURRENT_BINARY_DIR}/${HALMD_ARCHIVE_INTERMEDIATE_OUTPUT}
    --prefix=${HALMD_ARCHIVE_PREFIX}/
    ${HALMD_GIT_COMMIT_TAG}
  DEPENDS ${CMAKE_BINARY_DIR}/cmake/version.cmake
)

file(GLOB HALMD_ARCHIVE_LIBS RELATIVE ${CMAKE_SOURCE_DIR}/libs ${CMAKE_SOURCE_DIR}/libs/*)
list(SORT HALMD_ARCHIVE_LIBS)

foreach(lib ${HALMD_ARCHIVE_LIBS})
  set(HALMD_ARCHIVE_INTERMEDIATE_${lib} CMakeFiles/archive_${lib})

  add_custom_command(OUTPUT ${HALMD_ARCHIVE_INTERMEDIATE_${lib}}
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/libs/${lib}
    COMMAND git archive
      --format=tar
      --output=${CMAKE_CURRENT_BINARY_DIR}/${HALMD_ARCHIVE_INTERMEDIATE_${lib}}
      --prefix=${HALMD_ARCHIVE_PREFIX}/libs/${lib}/
      HEAD
    DEPENDS ${CMAKE_BINARY_DIR}/cmake/version.cmake
  )

  set(HALMD_ARCHIVE_INTERMEDIATE_INPUT ${HALMD_ARCHIVE_INTERMEDIATE_OUTPUT})
  set(HALMD_ARCHIVE_INTERMEDIATE_OUTPUT ${HALMD_ARCHIVE_INTERMEDIATE_INPUT}+${lib})

  add_custom_command(OUTPUT ${HALMD_ARCHIVE_INTERMEDIATE_OUTPUT}
    COMMAND cmake -E rename
      ${HALMD_ARCHIVE_INTERMEDIATE_INPUT}
      ${HALMD_ARCHIVE_INTERMEDIATE_OUTPUT}
    COMMAND tar
      --concatenate
      --file=${HALMD_ARCHIVE_INTERMEDIATE_OUTPUT}
      ${HALMD_ARCHIVE_INTERMEDIATE_${lib}}
    COMMAND cmake -E remove
      ${HALMD_ARCHIVE_INTERMEDIATE_${lib}}
    DEPENDS
      ${HALMD_ARCHIVE_INTERMEDIATE_INPUT}
      ${HALMD_ARCHIVE_INTERMEDIATE_${lib}}
  )
endforeach()

set(HALMD_ARCHIVE_INTERMEDIATE_INPUT ${HALMD_ARCHIVE_INTERMEDIATE_OUTPUT})
set(HALMD_ARCHIVE_INTERMEDIATE_OUTPUT ${HALMD_ARCHIVE_INTERMEDIATE_INPUT}+doc)

add_custom_command(OUTPUT ${HALMD_ARCHIVE_INTERMEDIATE_OUTPUT}
  COMMAND cmake -E rename
    ${HALMD_ARCHIVE_INTERMEDIATE_INPUT}
    ${HALMD_ARCHIVE_INTERMEDIATE_OUTPUT}
  COMMAND tar
    --append
    --file=${HALMD_ARCHIVE_INTERMEDIATE_OUTPUT}
    --owner=root
    --group=root
    --exclude=.doctrees
    --exclude=.buildinfo
    --transform="s,^doc/html/_sources,doc,"
    --transform="s,^doc/pdf,doc,"
    --transform="s,^,${HALMD_ARCHIVE_PREFIX}/,"
    "cmake/version.cmake"
    "doc/html"
    "doc/man/${PROGRAM_NAME}.1"
    "doc/pdf/${PROGRAM_NAME}.pdf"
  DEPENDS
    ${HALMD_ARCHIVE_INTERMEDIATE_INPUT}
    halmd_doc_html
    halmd_doc_man
    halmd_doc_pdf
)

add_custom_command(OUTPUT ${HALMD_ARCHIVE_OUTPUT}
  COMMAND cmake -E rename
    ${HALMD_ARCHIVE_INTERMEDIATE_OUTPUT}
    ${HALMD_ARCHIVE_PREFIX}.tar
  COMMAND bzip2
    --force ${HALMD_ARCHIVE_PREFIX}.tar
  DEPENDS
    ${HALMD_ARCHIVE_INTERMEDIATE_OUTPUT}
)

add_custom_target(halmd_archive
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  DEPENDS "${HALMD_ARCHIVE_PREFIX}.tar.bz2"
)
