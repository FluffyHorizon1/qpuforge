# Source stamp of the review and study programs, taken at build time (examples/git_stamp.hpp).
#
#   cmake -DSOURCE_DIR=<source tree> -DOUTPUT=<git_stamp.cpp> [-DGIT_EXECUTABLE=<git>] -P git_stamp.cmake
#
# Writes OUTPUT, a C++ source that defines review::git_stamp_commit and review::git_stamp_dirty:
#   commit  the full hash of HEAD (git rev-parse --verify HEAD), if SOURCE_DIR is the top level of a git work tree;
#   dirty   true if tracked files differ from that commit (git status --porcelain --untracked-files=no prints
#           something). Untracked files do not count.
# Without git, outside a git checkout (a source tarball, also one unpacked inside another repository) or when git
# fails, the stamp is "unknown" / false, and the script still succeeds.
#
# OUTPUT is rewritten only when the stamp changes, and one line is printed when it does. The target that runs this
# script on every build (qpuforge_git_stamp_update in CMakeLists.txt) therefore causes a recompile and relinks only
# then. GIT_EXECUTABLE is the git that CMake found when it configured the build directory; if that file is missing,
# git is looked up on PATH now. git runs with its optional locks off, so the build never writes to the repository.
cmake_minimum_required(VERSION 3.16)
foreach(v SOURCE_DIR OUTPUT)
  if(NOT DEFINED ${v} OR "${${v}}" STREQUAL "")
    message(FATAL_ERROR "git_stamp.cmake: -D${v}=... is required")
  endif()
endforeach()

set(git "${GIT_EXECUTABLE}")
if(NOT git OR NOT EXISTS "${git}")
  find_program(GIT_STAMP_GIT NAMES git PATHS ENV PATH NO_DEFAULT_PATH)
  set(git "${GIT_STAMP_GIT}")
endif()

set(commit "unknown")
set(dirty "false")
set(rc "")
set(why "git not found")   # then git's own message, for the line printed below
if(git)
  set(ENV{GIT_OPTIONAL_LOCKS} 0)
  execute_process(COMMAND "${git}" rev-parse --show-prefix WORKING_DIRECTORY "${SOURCE_DIR}"
                  RESULT_VARIABLE rc OUTPUT_VARIABLE prefix ERROR_VARIABLE why OUTPUT_STRIP_TRAILING_WHITESPACE)
  # An empty prefix: SOURCE_DIR is the top level. In a directory below it, git answers for the enclosing repository.
  if("${rc}" STREQUAL "0" AND "${prefix}" STREQUAL "")
    execute_process(COMMAND "${git}" rev-parse --verify HEAD WORKING_DIRECTORY "${SOURCE_DIR}"
                    RESULT_VARIABLE rc OUTPUT_VARIABLE head ERROR_VARIABLE why OUTPUT_STRIP_TRAILING_WHITESPACE)
    string(LENGTH "${head}" length)
    if("${rc}" STREQUAL "0" AND head MATCHES "^[0-9a-f]+$" AND (length EQUAL 40 OR length EQUAL 64))
      execute_process(COMMAND "${git}" status --porcelain --untracked-files=no WORKING_DIRECTORY "${SOURCE_DIR}"
                      RESULT_VARIABLE rc OUTPUT_VARIABLE status ERROR_VARIABLE why OUTPUT_STRIP_TRAILING_WHITESPACE)
      if("${rc}" STREQUAL "0")
        set(commit "${head}")
        if(NOT "${status}" STREQUAL "")
          set(dirty "true")
        endif()
      endif()
    endif()
  endif()
endif()

set(text "// Written by cmake/git_stamp.cmake, which runs on every build and rewrites this file when the stamp changes.
// Do not edit.
#include \"git_stamp.hpp\"
namespace review {
const char* const git_stamp_commit = \"${commit}\";
const bool git_stamp_dirty = ${dirty};
}
")
set(old "")
if(EXISTS "${OUTPUT}")
  file(READ "${OUTPUT}" old)
endif()
if(NOT "${old}" STREQUAL "${text}")
  file(WRITE "${OUTPUT}" "${text}")
  if(NOT commit STREQUAL "unknown")
    if(dirty)
      message(STATUS "source stamp: ${commit}, tracked files modified")
    else()
      message(STATUS "source stamp: ${commit}")
    endif()
  elseif(EXISTS "${SOURCE_DIR}/.git")   # a checkout whose commit cannot be read: say why
    string(REGEX REPLACE "[\r\n].*$" "" why "${why}")
    if("${why}" STREQUAL "")
      set(why "git: ${rc}")
    endif()
    message(STATUS "source stamp: unknown (${why})")
  else()
    message(STATUS "source stamp: unknown (not a git checkout)")
  endif()
endif()
