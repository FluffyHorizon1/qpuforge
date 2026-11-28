# Stamp guard (not one of the validation checks): cmake/git_stamp.cmake must describe the tree it is given at the
# moment it runs. This script makes a scratch git repository in WORK_DIR/git_stamp_test, runs the stamp script after
# each step and reads the stamp back from the file it writes:
#   git hidden from the script (no GIT_EXECUTABLE, an empty PATH)      "unknown" / false
#   no repository                                                      "unknown" / false
#   a repository without a commit                                      "unknown" / false
#   first commit                                                       git rev-parse HEAD, clean
#   git named by GIT_EXECUTABLE alone, then found on PATH alone        the same
#   the same again                                                     the stamp file is not rewritten
#   a tracked file touched, not changed                                the same; .git/index is not rewritten
#   an untracked file added                                            the same commit, still clean
#   a directory inside the work tree (an unpacked tarball)             "unknown" / false
#   second commit                                                      the new HEAD, clean
#   a tracked file modified, then staged                               the same commit, dirty
#   the modification committed                                         the new HEAD, clean again
#   git hidden again, in the checkout                                  "unknown" / false
# No network, and nothing outside WORK_DIR/git_stamp_test is written. Without git only the first check runs, and
# the "git_stamp: skipped" line makes CTest report the test as skipped (SKIP_REGULAR_EXPRESSION in CMakeLists.txt).
# Whether the build runs the stamp script at the right time is checked by tests/check_git_stamp_wiring.cmake.
#   cmake -DSTAMP_SCRIPT=<cmake/git_stamp.cmake> -DWORK_DIR=<build tree> [-DGIT_EXECUTABLE=<git>]
#         -P check_git_stamp.cmake
cmake_minimum_required(VERSION 3.16)
foreach(v STAMP_SCRIPT WORK_DIR)
  if(NOT DEFINED ${v} OR "${${v}}" STREQUAL "")
    message(FATAL_ERROR "check_git_stamp.cmake: -D${v}=... is required")
  endif()
endforeach()

set(work "${WORK_DIR}/git_stamp_test")
include("${CMAKE_CURRENT_LIST_DIR}/git_scratch.cmake")   # the scratch directory, `git`, check_failed, give_up, git_in
set(repo "${work}/repo")
set(stamp_cpp "${work}/out/git_stamp.cpp")   # out/ does not exist yet: the stamp script creates it
file(MAKE_DIRECTORY "${repo}" "${work}/no-git")

# Runs the stamp script on <dir> and compares the stamp it leaves in the stamp file with <commit> / <dirty>.
# stamp_git is the GIT_EXECUTABLE given to the script, stamp_env its environment changes (for cmake -E env).
function(expect_stamp what dir commit dirty)
  execute_process(COMMAND "${CMAKE_COMMAND}" -E env ${stamp_env} "${CMAKE_COMMAND}" "-DSOURCE_DIR=${dir}"
                          "-DOUTPUT=${stamp_cpp}" "-DGIT_EXECUTABLE=${stamp_git}" -P "${STAMP_SCRIPT}"
                  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE out)
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
    check_failed("${what}: the stamp script failed (${rc}):" "${out}")
  elseif(NOT "${got}" STREQUAL "${commit} / ${dirty}")
    check_failed("${what}: stamp ${got}, expected ${commit} / ${dirty}. The stamp script printed:" "${out}")
  else()
    message(STATUS "ok: ${what}: ${got}")
  endif()
endfunction()

# Gives <file> a date in 2000, where there is a touch to do it with. <result> is set to whether it was done.
find_program(GIT_STAMP_TEST_TOUCH NAMES touch PATHS ENV PATH NO_DEFAULT_PATH)
function(give_old_date file result)
  set(rc "touch not found")
  if(GIT_STAMP_TEST_TOUCH)
    execute_process(COMMAND "${GIT_STAMP_TEST_TOUCH}" -t 200006150000 "${file}" RESULT_VARIABLE rc ERROR_QUIET)
  endif()
  if("${rc}" STREQUAL "0")
    set(${result} TRUE PARENT_SCOPE)
  else()
    set(${result} FALSE PARENT_SCOPE)
  endif()
endfunction()

set(stamp_git "")
set(stamp_env "PATH=${work}/no-git")
expect_stamp("git hidden, no repository" "${repo}" unknown false)
if(NOT git)
  skip_without_git("git_stamp: skipped, git not found (the checks that need a repository were not run)")
  return()
endif()
set(stamp_git "${git}")
set(stamp_env "")

# With the build tree inside the qpuforge checkout, this is also a directory inside another repository's work tree.
expect_stamp("no repository" "${repo}" unknown false)

git_in("${repo}" init -q)
file(WRITE "${repo}/tracked.txt" "one\n")
git_in("${repo}" add tracked.txt)
expect_stamp("repository without a commit" "${repo}" unknown false)

git_in("${repo}" commit -q -m first)
git_in("${repo}" rev-parse HEAD)
set(first "${git_out}")
expect_stamp("first commit" "${repo}" "${first}" false)

# Which git the script runs: GIT_EXECUTABLE if that file exists, also when PATH has no git (CMake was pointed at one
# when it configured); otherwise the git on PATH (git was installed, or moved, after the build directory was
# configured). Each of the three starts without a stamp file, so none can pass on what the last one wrote.
get_filename_component(git_directory "${git}" DIRECTORY)
file(REMOVE "${stamp_cpp}")
set(stamp_env "PATH=${work}/no-git")
expect_stamp("git named by GIT_EXECUTABLE, none on PATH" "${repo}" "${first}" false)
file(REMOVE "${stamp_cpp}")
set(stamp_git "")
set(stamp_env "PATH=${git_directory}")
expect_stamp("no GIT_EXECUTABLE, git on PATH" "${repo}" "${first}" false)
file(REMOVE "${stamp_cpp}")
set(stamp_git "${work}/no-git/git")
expect_stamp("GIT_EXECUTABLE names a missing file, git on PATH" "${repo}" "${first}" false)
set(stamp_git "${git}")
set(stamp_env "")

# Nothing changed: the stamp file must be left alone, or every build would recompile it and relink the programs.
# The file is given an old date; a rewrite would bring it back to today.
give_old_date("${stamp_cpp}" old)
if(old)
  expect_stamp("first commit again" "${repo}" "${first}" false)
  file(TIMESTAMP "${stamp_cpp}" year "%Y")
  if("${year}" STREQUAL "2000")
    message(STATUS "ok: the stamp file was not rewritten")
  else()
    check_failed("the stamp file was rewritten although the stamp had not changed (its year is now ${year})")
  endif()
else()
  message(STATUS "not checked: that an unchanged stamp file is not rewritten (no touch -t)")
endif()

# The script must not write to the repository it reads. A tracked file gets another date and the same content: git
# status then refreshes the index in memory, and writes it back to .git/index unless its optional locks are off
# (GIT_OPTIONAL_LOCKS, which git has had since 2.15).
git_in("${repo}" --version)
if(git_out MATCHES "^git version ([0-9]+\\.[0-9]+)" AND NOT CMAKE_MATCH_1 VERSION_LESS 2.15)
  file(WRITE "${repo}/tracked.txt" "one\n")
  give_old_date("${repo}/tracked.txt" old)
  file(SHA256 "${repo}/.git/index" index_before)
  expect_stamp("tracked file touched, not changed" "${repo}" "${first}" false)
  file(SHA256 "${repo}/.git/index" index_after)
  if("${index_after}" STREQUAL "${index_before}")
    message(STATUS "ok: .git/index was not rewritten")
  else()
    check_failed("the stamp script rewrote .git/index: git must run with its optional locks off")
  endif()
else()
  message(STATUS "not checked: that .git/index is not rewritten (${git_out})")
endif()

file(WRITE "${repo}/untracked.txt" "not added\n")
expect_stamp("untracked file added" "${repo}" "${first}" false)

file(MAKE_DIRECTORY "${repo}/unpacked")
expect_stamp("directory inside the work tree" "${repo}/unpacked" unknown false)

file(APPEND "${repo}/tracked.txt" "two\n")
git_in("${repo}" commit -q -a -m second)
git_in("${repo}" rev-parse HEAD)
set(second "${git_out}")
if("${second}" STREQUAL "${first}")
  give_up("the second commit has the hash of the first")
endif()
expect_stamp("second commit" "${repo}" "${second}" false)

file(APPEND "${repo}/tracked.txt" "three\n")
expect_stamp("tracked file modified" "${repo}" "${second}" true)
git_in("${repo}" add tracked.txt)
expect_stamp("modification staged" "${repo}" "${second}" true)

git_in("${repo}" commit -q -m third)
git_in("${repo}" rev-parse HEAD)
set(third "${git_out}")
expect_stamp("modification committed" "${repo}" "${third}" false)

set(stamp_git "")
set(stamp_env "PATH=${work}/no-git")
expect_stamp("git hidden, in the checkout" "${repo}" unknown false)

file(REMOVE_RECURSE "${work}")
