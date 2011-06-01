# - GitRepository
# Extract commit metadata from Git repository

#=============================================================================
# Copyright 2010-2011 Peter Colberg
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file COPYING-CMAKE-SCRIPTS for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================

# This macro defines
#  <PREFIX>_GIT_BRANCH
#  <PREFIX>_GIT_COMMIT_TAG
#  <PREFIX>_GIT_COMMITTER_DATE
#
macro(git_repository DIR PREFIX)

  # Find path to .git directory. If <DIR> is the top-level directory,
  # _GIT_DIR will contain the relative path ".git". If <DIR> is a
  # subdirectory, _GIT_DIR will contain an absolute path.
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" rev-parse --git-dir
    WORKING_DIRECTORY "${DIR}"
    RESULT_VARIABLE _GIT_STATUS
    OUTPUT_VARIABLE _GIT_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT _GIT_STATUS EQUAL 0)
    message(FATAL_ERROR "git_repository Failure finding repository")
  endif()
  if(NOT IS_ABSOLUTE "${_GIT_DIR}")
    set(_GIT_DIR "${DIR}/${_GIT_DIR}")
  endif()

  # Rerun CMake if .git/logs/HEAD is touched, which occurs if commit changes.
  configure_file(
    "${_GIT_DIR}/logs/HEAD"
    "CMakeFiles/GitRepository/${PREFIX}/logs/HEAD"
    COPYONLY
  )

  # A detached HEAD will result in a return value equal to 1.
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" symbolic-ref -q HEAD
    WORKING_DIRECTORY "${DIR}"
    RESULT_VARIABLE _GIT_STATUS
    OUTPUT_VARIABLE _GIT_SYMBOLIC_REF
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT _GIT_STATUS EQUAL 0 AND NOT _GIT_STATUS EQUAL 1)
    message(FATAL_ERROR "git_repository Failure executing git symbolic-ref")
  endif()

  if(_GIT_SYMBOLIC_REF MATCHES "^refs/heads/")
    string(REGEX REPLACE "^refs/heads/" "" ${PREFIX}_GIT_BRANCH "${_GIT_SYMBOLIC_REF}")
  endif()

  # Parse commit tag
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" describe --always HEAD
    WORKING_DIRECTORY "${DIR}"
    RESULT_VARIABLE _GIT_STATUS
    OUTPUT_VARIABLE ${PREFIX}_GIT_COMMIT_TAG
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT _GIT_STATUS EQUAL 0)
    message(FATAL_ERROR "git_repository Failure executing git describe")
  endif()

  # Parse committer date
  execute_process(
    COMMAND "${GIT_EXECUTABLE}" log -1 --pretty=format:%cD HEAD
    WORKING_DIRECTORY "${DIR}"
    RESULT_VARIABLE _GIT_STATUS
    OUTPUT_VARIABLE ${PREFIX}_GIT_COMMITTER_DATE
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(NOT _GIT_STATUS EQUAL 0)
    message(FATAL_ERROR "git_repository Failure executing git log")
  endif()

endmacro()
