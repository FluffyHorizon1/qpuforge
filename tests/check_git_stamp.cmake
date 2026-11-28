# Stamp guard (not one of the validation checks): cmake/git_stamp.cmake must describe the tree it is given at the
# moment it runs. This script makes a scratch git repository in WORK_DIR/git_stamp_test, runs the stamp script after
# each step and reads the stamp back from the file it writes:
#   git hidden from the script (no GIT_EXECUTABLE, an empty PATH)      "unknown" / false
#   no repository                                                      "unknown" / false
#   a repository without a commit                                      "unknown" / false
#   first commit                                                       git rev-parse HEAD, clean
#   the same again                                                     the stamp file is not rewritten
#   an untracked file added                                            the same commit, still clean
#   a directory inside the work tree (an unpacked tarball)             "unknown" / false
#   second commit                                                      the new HEAD, clean
#   a tracked file modified, then staged                               the same commit, dirty
#   the modification committed                                         the new HEAD, clean again
#   git hidden again, in the checkout                                  "unknown" / false
# No network, and nothing outside WORK_DIR/git_stamp_test is written. Without git only the first check runs, and
# the "git_stamp: skipped" line makes CTest report the test as skipped (SKIP_REGULAR_EXPRESSION in CMakeLists.txt).
#   cmake -DSTAMP_SCRIPT=<cmake/git_stamp.cmake> -DWORK_DIR=<build tree> [-DGIT_EXECUTABLE=<git>]
#         -P check_git_stamp.cmake
cmake_minimum_required(VERSION 3.16)
foreach(v STAMP_SCRIPT WORK_DIR)
  if(NOT DEFINED ${v} OR "${${v}}" STREQUAL "")
    message(FATAL_ERROR "check_git_stamp.cmake: -D${v}=... is required")
  endif()
endforeach()

set(work "${WORK_DIR}/git_stamp_test")
set(repo "${work}/repo")
set(stamp_cpp "${work}/out/git_stamp.cpp")   # out/ does not exist yet: the stamp script creates it
file(REMOVE_RECURSE "${work}")
file(MAKE_DIRECTORY "${repo}" "${work}/home" "${work}/no-git")

# The scratch repository must not pick up the user's git configuration (identity, commit signing, hooks, templates),
# and git must not be pointed at another repository by the environment.
foreach(v GIT_DIR GIT_WORK_TREE GIT_INDEX_FILE GIT_OBJECT_DIRECTORY GIT_ALTERNATE_OBJECT_DIRECTORIES GIT_COMMON_DIR
          GIT_NAMESPACE GIT_CEILING_DIRECTORIES GIT_TEMPLATE_DIR GIT_CONFIG GIT_CONFIG_COUNT GIT_CONFIG_PARAMETERS)
  unset(ENV{${v}})
endforeach()
file(WRITE "${work}/home/.gitconfig" "")
set(ENV{HOME} "${work}/home")
set(ENV{XDG_CONFIG_HOME} "${work}/home/.config")
set(ENV{GIT_CONFIG_GLOBAL} "${work}/home/.gitconfig")
set(ENV{GIT_CONFIG_NOSYSTEM} 1)
foreach(who AUTHOR COMMITTER)
  set(ENV{GIT_${who}_NAME} "qpuforge stamp test")
  set(ENV{GIT_${who}_EMAIL} "stamp-test@example.invalid")
endforeach()

# the same choice of git as cmake/git_stamp.cmake
set(git "${GIT_EXECUTABLE}")
if(NOT git OR NOT EXISTS "${git}")
  find_program(GIT_STAMP_TEST_GIT NAMES git PATHS ENV PATH NO_DEFAULT_PATH)
  set(git "${GIT_STAMP_TEST_GIT}")
endif()

# A failure of the scaffolding (not of a check): remove the scratch repository, then stop.
function(give_up why)
  file(REMOVE_RECURSE "${work}")
  message(FATAL_ERROR "check_git_stamp.cmake: ${why}")
endfunction()

# git <arguments> in the scratch repository; its output is left in git_out.
function(git_in_repo)
  execute_process(COMMAND "${git}" ${ARGN} WORKING_DIRECTORY "${repo}"
                  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT "${rc}" STREQUAL "0")
    give_up("git ${ARGN} failed (${rc}):\n${out}\n${err}")
  endif()
  set(git_out "${out}" PARENT_SCOPE)
endfunction()

# Runs the stamp script on <dir> and compares the stamp it leaves in the stamp file with <commit> / <dirty>.
# stamp_git is the GIT_EXECUTABLE given to the script, stamp_env its environment changes (for cmake -E env).
function(expect_stamp what dir commit dirty)
  execute_process(COMMAND "${CMAKE_COMMAND}" -E env ${stamp_env} "${CMAKE_COMMAND}" "-DSOURCE_DIR=${dir}"
                          "-DOUTPUT=${stamp_cpp}" "-DGIT_EXECUTABLE=${stamp_git}" -P "${STAMP_SCRIPT}"
                  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
  set(text "")
  if(EXISTS "${stamp_cpp}")
    file(READ "${stamp_cpp}" text)
  endif()
  set(got "no stamp")
  if(text MATCHES "git_stamp_commit = \"([^\"]*)\";")
    set(got "${CMAKE_MATCH_1}")
    if(text MATCHES "git_stamp_dirty = (true|false);")
      set(got "${got} / ${CMAKE_MATCH_1}")
    endif()
  endif()
  if(NOT "${rc}" STREQUAL "0")
    message(SEND_ERROR "${what}: the stamp script failed (${rc}):\n${out}\n${err}")
  elseif(NOT "${got}" STREQUAL "${commit} / ${dirty}")
    message(SEND_ERROR "${what}: stamp ${got}, expected ${commit} / ${dirty}\n${out}\n${err}")
  else()
    message(STATUS "ok: ${what}: ${got}")
  endif()
endfunction()

set(stamp_git "")
set(stamp_env "PATH=${work}/no-git")
expect_stamp("git hidden, no repository" "${repo}" unknown false)
if(NOT git)
  file(REMOVE_RECURSE "${work}")
  message(STATUS "git_stamp: skipped, git not found (the checks that need a repository were not run)")
  return()
endif()
set(stamp_git "${git}")
set(stamp_env "")

# With the build tree inside the qpuforge checkout, this is also a directory inside another repository's work tree.
expect_stamp("no repository" "${repo}" unknown false)

git_in_repo(init -q)
file(WRITE "${repo}/tracked.txt" "one\n")
git_in_repo(add tracked.txt)
expect_stamp("repository without a commit" "${repo}" unknown false)

git_in_repo(commit -q -m first)
git_in_repo(rev-parse HEAD)
set(first "${git_out}")
expect_stamp("first commit" "${repo}" "${first}" false)

# Nothing changed: the stamp file must be left alone, or every build would recompile it and relink the programs.
# The file is given an old date; a rewrite would bring it back to today.
find_program(GIT_STAMP_TEST_TOUCH NAMES touch PATHS ENV PATH NO_DEFAULT_PATH)
set(rc "touch not found")
if(GIT_STAMP_TEST_TOUCH)
  execute_process(COMMAND "${GIT_STAMP_TEST_TOUCH}" -t 200006150000 "${stamp_cpp}" RESULT_VARIABLE rc ERROR_QUIET)
endif()
if("${rc}" STREQUAL "0")
  expect_stamp("first commit again" "${repo}" "${first}" false)
  file(TIMESTAMP "${stamp_cpp}" year "%Y")
  if("${year}" STREQUAL "2000")
    message(STATUS "ok: the stamp file was not rewritten")
  else()
    message(SEND_ERROR "the stamp file was rewritten although the stamp had not changed (its year is now ${year})")
  endif()
else()
  message(STATUS "not checked: that an unchanged stamp file is not rewritten (touch -t: ${rc})")
endif()

file(WRITE "${repo}/untracked.txt" "not added\n")
expect_stamp("untracked file added" "${repo}" "${first}" false)

file(MAKE_DIRECTORY "${repo}/unpacked")
expect_stamp("directory inside the work tree" "${repo}/unpacked" unknown false)

file(APPEND "${repo}/tracked.txt" "two\n")
git_in_repo(commit -q -a -m second)
git_in_repo(rev-parse HEAD)
set(second "${git_out}")
if("${second}" STREQUAL "${first}")
  give_up("the second commit has the hash of the first")
endif()
expect_stamp("second commit" "${repo}" "${second}" false)

file(APPEND "${repo}/tracked.txt" "three\n")
expect_stamp("tracked file modified" "${repo}" "${second}" true)
git_in_repo(add tracked.txt)
expect_stamp("modification staged" "${repo}" "${second}" true)

git_in_repo(commit -q -m third)
git_in_repo(rev-parse HEAD)
set(third "${git_out}")
expect_stamp("modification committed" "${repo}" "${third}" false)

set(stamp_git "")
set(stamp_env "PATH=${work}/no-git")
expect_stamp("git hidden, in the checkout" "${repo}" unknown false)

file(REMOVE_RECURSE "${work}")
