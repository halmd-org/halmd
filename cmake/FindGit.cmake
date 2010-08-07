# - Extract information from a Git repository
# The module defines the following variables:
#  GIT_EXECUTABLE - path to git executable
#  GIT_VERSION - git version
#  GIT_FOUND - true if git was found
# If the git executable is found the macro
#  GIT_REPOSITORY_INFO(<dir> <var-prefix>)
# is defined to extract information of a git repository at
# a given location. The macro defines the following variables:
#  <var-prefix>_GIT_COMMIT - commit hash of current branch tip
#  <var-prefix>_GIT_REFS - ref names
#  <var-prefix>_GIT_TAGS - tags
#  <var-prefix>_GIT_HEADS - local branches
#  <var-prefix>_GIT_REMOTES - remote branches
#  <var-prefix>_GIT_AUTHOR - author name and email
#  <var-prefix>_GIT_AUTHOR_DATE - author date
#  <var-prefix>_GIT_COMMITTER - committer name and email
#  <var-prefix>_GIT_COMMITTER_DATE - committer date
#  <var-prefix>_GIT_SUBJECT - subject
#  <var-prefix>_GIT_VERSION - most recent tag
# Example usage:
#  FIND_PACKAGE(Git)
#  if(GIT_FOUND)
#    GIT_REPOSITORY_INFO(${PROJECT_SOURCE_DIR} Project)
#    MESSAGE("Current revision is ${Project_GIT_COMMIT}")
#  endif(GIT_FOUND)

# Copyright (C) 2008  Peter Colberg
#
# This file was derived from FindSubversion.cmake shipped with CMake 2.4.7.
#
# Copyright (c) 2006, Tristan Carel
# All rights reserved.
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the University of California, Berkeley nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE REGENTS AND CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN if ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


find_program(GIT_EXECUTABLE git
  DOC "Git version control system executable")
mark_as_advanced(GIT_EXECUTABLE)

if(GIT_EXECUTABLE)
  set(GIT_FOUND TRUE)

  macro(GIT_COMMAND dir command)
    execute_process(COMMAND ${GIT_EXECUTABLE} ${command} ${ARGN}
      WORKING_DIRECTORY ${dir}
      OUTPUT_VARIABLE GIT_${command}_OUTPUT
      ERROR_VARIABLE GIT_${command}_ERROR
      RESULT_VARIABLE GIT_${command}_RESULT
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    if(NOT ${GIT_${command}_RESULT} EQUAL 0)
      set(cmdline "${GIT_EXECUTABLE} ${command}")
      foreach(arg ${ARGN})
	set(cmdline "${cmdline} ${arg}")
      endforeach(arg ${ARGN})
      message(SEND_ERROR "Command \"${cmdline}\" failed with output:\n${GIT_${command}_ERROR}")

      set(GIT_${command}_OUTPUT)
    endif(NOT ${GIT_${command}_RESULT} EQUAL 0)

  endmacro(GIT_COMMAND dir command)

  # git version
  GIT_COMMAND(${PROJECT_SOURCE_DIR} version)
  string(REGEX REPLACE "^git version ([.0-9]+).*" "\\1" GIT_VERSION "${GIT_version_OUTPUT}")

  macro(GIT_REPOSITORY_INFO dir prefix)
    if(IS_DIRECTORY "${dir}")
      # generate commit log for current branch tip
      GIT_COMMAND(${dir} log -1 --decorate --pretty=fuller)

      # commit hash
      string(REGEX REPLACE "^commit ([0-9a-f]+).*" "\\1" ${prefix}_GIT_COMMIT "${GIT_log_OUTPUT}")

      # ref names
      if(GIT_log_OUTPUT MATCHES "^commit [0-9a-f]+ \\([^\n]+\\).*")
	string(REGEX REPLACE "^commit [0-9a-f]+ \\(([^\n]+)\\).*" "\\1" ${prefix}_GIT_REFS "${GIT_log_OUTPUT}")
	string(REGEX REPLACE ", " ";" ${prefix}_GIT_REFS "${${prefix}_GIT_REFS}")
      endif(GIT_log_OUTPUT MATCHES "^commit [0-9a-f]+ \\([^\n]+\\).*")

      # author name and email
      if(GIT_log_OUTPUT MATCHES ".*\nAuthor:([^\n]+).*")
	string(REGEX REPLACE ".*\nAuthor:([^\n]+).*" "\\1" ${prefix}_GIT_AUTHOR "${GIT_log_OUTPUT}")
	string(STRIP "${${prefix}_GIT_AUTHOR}" ${prefix}_GIT_AUTHOR)
      endif(GIT_log_OUTPUT MATCHES ".*\nAuthor:([^\n]+).*")

      # author date
      if(GIT_log_OUTPUT MATCHES ".*\nAuthorDate:([^\n]+).*")
	string(REGEX REPLACE ".*\nAuthorDate:([^\n]+).*" "\\1" ${prefix}_GIT_AUTHOR_DATE "${GIT_log_OUTPUT}")
	string(STRIP "${${prefix}_GIT_AUTHOR_DATE}" ${prefix}_GIT_AUTHOR_DATE)
      endif(GIT_log_OUTPUT MATCHES ".*\nAuthorDate:([^\n]+).*")

      # committer name and email
      if(GIT_log_OUTPUT MATCHES ".*\nCommit:([^\n]+).*")
	string(REGEX REPLACE ".*\nCommit:([^\n]+).*" "\\1" ${prefix}_GIT_COMMITTER "${GIT_log_OUTPUT}")
	string(STRIP "${${prefix}_GIT_COMMITTER}" ${prefix}_GIT_COMMITTER)
      endif(GIT_log_OUTPUT MATCHES ".*\nCommit:([^\n]+).*")

      # committer date
      if(GIT_log_OUTPUT MATCHES ".*\nCommitDate:([^\n]+).*")
	string(REGEX REPLACE ".*\nCommitDate:([^\n]+).*" "\\1" ${prefix}_GIT_COMMITTER_DATE "${GIT_log_OUTPUT}")
	string(STRIP "${${prefix}_GIT_COMMITTER_DATE}" ${prefix}_GIT_COMMITTER_DATE)
      endif(GIT_log_OUTPUT MATCHES ".*\nCommitDate:([^\n]+).*")

      # subject
      if(GIT_log_OUTPUT MATCHES ".*\n\n([^\n]+).*")
	string(REGEX REPLACE ".*\n\n([^\n]+).*" "\\1" ${prefix}_GIT_SUBJECT "${GIT_log_OUTPUT}")
	string(STRIP "${${prefix}_GIT_SUBJECT}" ${prefix}_GIT_SUBJECT)
      endif(GIT_log_OUTPUT MATCHES ".*\n\n([^\n]+).*")

      foreach(ref ${${prefix}_GIT_REFS})
	# tags
	if(ref MATCHES "^refs/tags/.*$")
	  string(REGEX REPLACE "^refs/tags/(.*)$" "\\1" tag "${ref}")
	  list(APPEND ${prefix}_GIT_TAGS "${tag}")
	endif(ref MATCHES "^refs/tags/.*$")

	# local branches
	if(ref MATCHES "^refs/heads/.*$")
	  string(REGEX REPLACE "^refs/heads/(.*)$" "\\1" head "${ref}")
	  list(APPEND ${prefix}_GIT_HEADS "${head}")
	endif(ref MATCHES "^refs/heads/.*$")

	# remote branches
	if(ref MATCHES "^refs/remotes/.*$")
	  string(REGEX REPLACE "^refs/remotes/(.*)$" "\\1" remote "${ref}")
	  list(APPEND ${prefix}_GIT_REMOTES "${remote}")
	endif(ref MATCHES "^refs/remotes/.*$")
      endforeach(ref)

      # show most recent tag that is reachable from a commit
      GIT_COMMAND(${dir} describe --long --always)

      if(GIT_describe_OUTPUT)
	string(STRIP "${GIT_describe_OUTPUT}" ${prefix}_GIT_VERSION)
      endif(GIT_describe_OUTPUT)

    else(IS_DIRECTORY "${dir}")
      message(SEND_ERROR "Invalid GIT_REPOSITORY_INFO directory \"${dir}\"")
    endif(IS_DIRECTORY "${dir}")

  endmacro(GIT_REPOSITORY_INFO dir prefix)
endif(GIT_EXECUTABLE)

if(NOT GIT_FOUND)
  if(NOT GIT_FIND_QUIETLY)
    MESSAGE(STATUS "Git executable was not found.")
  else(NOT GIT_FIND_QUIETLY)
    if(GIT_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Git executable was not found.")
    endif(GIT_FIND_REQUIRED)
  endif(NOT GIT_FIND_QUIETLY)
endif(NOT GIT_FOUND)
