# Wiring guard (not one of the validation checks): the build must take the source stamp each time it BUILDS, and every
# review and study program must carry it. tests/check_git_stamp.cmake checks the script cmake/git_stamp.cmake; this
# script checks CMakeLists.txt, which decides when that script runs and which programs are linked with its result.
#
# It makes a copy of the build description in WORK_DIR/git_stamp_wiring_test: CMakeLists.txt, cmake/ and
# examples/git_stamp.hpp as they are, and a stand-in for every program source. The stand-in of a program that
# includes git_stamp.hpp prints the stamp the program was built with, so nothing of qpuforge itself is compiled and a
# build takes about a second. The copy is configured ONCE, with the generator and the compiler of this build. Its
# programs are then built and run after each step:
#   no repository (a source tarball)             "unknown" / false
#   the copy made a repository, first commit     git rev-parse HEAD, clean
#   a second commit that changes no source       the new HEAD, clean; only the stamp is compiled
#   nothing changed                              the same; nothing is compiled or linked
#   a tracked file modified                      the same commit, dirty
# A stamp taken when CMake configures (the old behaviour) fails at the first commit. A program that CMakeLists.txt
# does not link with the stamp library prints "unknown" and fails there too. The builds must also run the git that
# the copy was configured with (GIT_EXECUTABLE has to reach cmake/git_stamp.cmake).
# No network, and nothing outside WORK_DIR/git_stamp_wiring_test is written. Without git only the first build runs,
# and the "git_stamp_wiring: skipped" line makes CTest report the test as skipped (SKIP_REGULAR_EXPRESSION in
# CMakeLists.txt).
#   cmake -DSOURCE_DIR=<source tree> -DWORK_DIR=<build tree> -DGENERATOR=<generator> -DCXX_COMPILER=<compiler>
#         [-DMAKE_PROGRAM=<make or ninja>] [-DGIT_EXECUTABLE=<git>] -P check_git_stamp_wiring.cmake
cmake_minimum_required(VERSION 3.16)
foreach(v SOURCE_DIR WORK_DIR GENERATOR CXX_COMPILER)
  if(NOT DEFINED ${v} OR "${${v}}" STREQUAL "")
    message(FATAL_ERROR "check_git_stamp_wiring.cmake: -D${v}=... is required")
  endif()
endforeach()

set(work "${WORK_DIR}/git_stamp_wiring_test")
include("${CMAKE_CURRENT_LIST_DIR}/git_scratch.cmake")   # the scratch directory, `git`, check_failed, give_up, git_in
set(src "${work}/src")
set(bld "${work}/build")

# Does <file> include git_stamp.hpp, itself or through the headers it includes from its own directory?
function(includes_git_stamp file result)
  get_filename_component(dir "${file}" DIRECTORY)
  set(todo "${file}")
  set(seen "")
  set(found FALSE)
  while(NOT "${todo}" STREQUAL "" AND NOT found)
    list(GET todo 0 current)
    list(REMOVE_AT todo 0)
    list(APPEND seen "${current}")
    file(STRINGS "${current}" lines REGEX "^[ \t]*#[ \t]*include[ \t]*\"[^\"]+\"")
    foreach(line IN LISTS lines)
      string(REGEX REPLACE "^[^\"]*\"([^\"]+)\".*$" "\\1" header "${line}")
      if(header STREQUAL "git_stamp.hpp")
        set(found TRUE)
      elseif(EXISTS "${dir}/${header}" AND NOT "${dir}/${header}" IN_LIST seen)
        list(APPEND todo "${dir}/${header}")
      endif()
    endforeach()
  endwhile()
  set(${result} ${found} PARENT_SCOPE)
endfunction()

# The copy. The programs that must carry the stamp are those that write a provenance block: every examples/<name>.cpp
# that includes git_stamp.hpp (through review_io.hpp), which CMakeLists.txt builds as the target <name>.
file(COPY "${SOURCE_DIR}/CMakeLists.txt" "${SOURCE_DIR}/cmake" DESTINATION "${src}" NO_SOURCE_PERMISSIONS)
file(COPY "${SOURCE_DIR}/examples/git_stamp.hpp" DESTINATION "${src}/examples" NO_SOURCE_PERMISSIONS)
file(GLOB sources RELATIVE "${SOURCE_DIR}"
     "${SOURCE_DIR}/src/*.cpp" "${SOURCE_DIR}/tests/*.cpp" "${SOURCE_DIR}/examples/*.cpp")
set(programs "")
set(control "")   # a program that has nothing to do with the stamp
foreach(source IN LISTS sources)
  set(stamped FALSE)
  if(source MATCHES "^examples/([^/]+)\\.cpp$")
    set(name "${CMAKE_MATCH_1}")
    includes_git_stamp("${SOURCE_DIR}/${source}" stamped)
    if(NOT stamped)
      set(control "${name}")
    endif()
  endif()
  if(stamped)
    list(APPEND programs "${name}")
    file(WRITE "${src}/${source}" "// Stand-in for ${source}, written by tests/check_git_stamp_wiring.cmake.
#include \"git_stamp.hpp\"
#include <cstdio>
int main() {
    std::printf(\"%s %s\\n\", review::git_stamp_commit, review::git_stamp_dirty ? \"true\" : \"false\");
    return 0;
}
")
  else()
    file(WRITE "${src}/${source}" "// Stand-in for ${source}, written by tests/check_git_stamp_wiring.cmake.
int main() { return 0; }
")
  endif()
endforeach()
if("${programs}" STREQUAL "")
  give_up("no examples/*.cpp in ${SOURCE_DIR} includes git_stamp.hpp: there is nothing to check")
endif()
list(LENGTH programs count)
string(REPLACE ";" " " names "${programs}")
message(STATUS "programs that include git_stamp.hpp: ${names}")
file(WRITE "${src}/notes.txt" "one\n")   # the tracked file that the later steps change; no target reads it

# The git of the copy. CMakeLists.txt must hand the git that CMake has when it configures (GIT_EXECUTABLE) on to
# cmake/git_stamp.cmake, which looks for git on PATH only if that file is gone. To tell which git the builds run, the
# copy is configured with a git of its own where a shell script can be one: it notes each call in ${work}/git-calls
# and then runs the git of this test.
set(git_of_the_copy "${git}")
set(git_calls "")
if(git AND CMAKE_HOST_UNIX AND NOT "${work}${git}" MATCHES "'")
  file(WRITE "${work}/script/git" "#!/bin/sh\necho called >> '${work}/git-calls'\nexec '${git}' \"$@\"\n")
  file(COPY "${work}/script/git" DESTINATION "${work}/noting" FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE)
  execute_process(COMMAND "${work}/noting/git" --version RESULT_VARIABLE rc OUTPUT_QUIET ERROR_QUIET)
  if("${rc}" STREQUAL "0")
    set(git_of_the_copy "${work}/noting/git")
    set(git_calls "${work}/git-calls")
  endif()
endif()

set(configure "${CMAKE_COMMAND}" -S "${src}" -B "${bld}" -G "${GENERATOR}" "-DCMAKE_CXX_COMPILER=${CXX_COMPILER}")
if(MAKE_PROGRAM)
  list(APPEND configure "-DCMAKE_MAKE_PROGRAM=${MAKE_PROGRAM}")
endif()
if(git)
  list(APPEND configure "-DGIT_EXECUTABLE=${git_of_the_copy}")
endif()
execute_process(COMMAND ${configure} RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE out)
if(NOT "${rc}" STREQUAL "0")
  give_up("CMake could not configure the copy of the build description (${rc}). If CMakeLists.txt now names a source outside src/, tests/ and examples/, this script has to write a stand-in for it." "${out}")
endif()
file(REMOVE "${work}/git-calls")   # CMake asked that git for its version; only the calls of the builds count

# make and ninja compare file times, so what the next build writes must be newer than what the last one made. This
# waits until the file system gives a later time: at once on most, a second or two where times are whole seconds.
function(wait_for_a_later_file_time)
  file(WRITE "${work}/built" "")
  foreach(i RANGE 1 300)
    file(WRITE "${work}/later" "${i}")
    if(NOT "${work}/built" IS_NEWER_THAN "${work}/later")   # IS_NEWER_THAN is also true for equal times
      return()
    endif()
    execute_process(COMMAND "${CMAKE_COMMAND}" -E sleep 0.01)
  endforeach()
endfunction()

# Builds the programs in the copy, which is never configured again, runs each and compares the stamp it prints with
# <commit> / <dirty>. Leaves what the build printed in build_out, and whether the check passed in build_ok.
function(build_and_expect what commit dirty)
  set(build_ok FALSE PARENT_SCOPE)
  execute_process(COMMAND "${CMAKE_COMMAND}" --build "${bld}" --parallel 2 --target ${programs}
                  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE out)
  set(build_out "${out}" PARENT_SCOPE)
  if(NOT "${rc}" STREQUAL "0")
    check_failed("${what}: the build failed (${rc}):" "${out}")
    return()
  endif()
  set(wrong "")
  foreach(name IN LISTS programs)
    find_program(GIT_STAMP_WIRING_${name} NAMES "${name}" PATHS "${bld}"
                 PATH_SUFFIXES Debug Release RelWithDebInfo MinSizeRel NO_DEFAULT_PATH)
    execute_process(COMMAND "${GIT_STAMP_WIRING_${name}}"
                    RESULT_VARIABLE rc OUTPUT_VARIABLE printed ERROR_VARIABLE printed OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT "${rc}" STREQUAL "0")
      string(APPEND wrong "${name} could not be run (${rc}) ${printed}\n")
    elseif(NOT "${printed}" STREQUAL "${commit} ${dirty}")
      string(APPEND wrong "${name} printed: ${printed}\n")
    endif()
  endforeach()
  if(NOT "${wrong}" STREQUAL "")
    check_failed("${what}: every program must print the stamp \"${commit} ${dirty}\", but" "${wrong}The build printed:\n${out}")
  else()
    message(STATUS "ok: ${what}: ${count} programs with the stamp ${commit} / ${dirty}")
    set(build_ok TRUE PARENT_SCOPE)
  endif()
  wait_for_a_later_file_time()
endfunction()

build_and_expect("no repository" unknown false)
# What the later builds compile and link is checked too, if it can be. The generator must say it, as Ninja and the
# Makefile generators do. And a source that has not changed must be left alone in this build tree, whatever the stamp
# does: to see that, a program that has nothing to do with the stamp is built twice. (Ninja 1.11.1 with gcc 13
# compiles every source again on every build if a path has an apostrophe in it.)
set(cost_not_checked "")   # the reason, if it cannot be
if(NOT build_out MATCHES "Building CXX object" OR NOT build_out MATCHES "Linking CXX executable")
  set(cost_not_checked "${GENERATOR} does not say what it compiles and links as Ninja and make do")
elseif(NOT "${control}" STREQUAL "")
  foreach(build first second)
    execute_process(COMMAND "${CMAKE_COMMAND}" --build "${bld}" --target "${control}"
                    RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE out)
  endforeach()
  if("${rc}" STREQUAL "0" AND out MATCHES "Building CXX object")
    set(cost_not_checked "every build here compiles ${control} again, a source that has not changed")
  endif()
endif()
if(NOT git)
  skip_without_git("git_stamp_wiring: skipped, git not found (the builds that need a repository were not run)")
  return()
endif()
if("${git_calls}" STREQUAL "")
  message(STATUS "not checked: which git the build runs (no shell script could stand in for git here)")
elseif(EXISTS "${git_calls}")
  message(STATUS "ok: the build ran the git that the copy was configured with")
else()
  check_failed("the build did not run the git that the copy was configured with: CMakeLists.txt must pass GIT_EXECUTABLE on to cmake/git_stamp.cmake")
endif()

git_in("${src}" init -q)
git_in("${src}" add -A)
git_in("${src}" commit -q -m first)
git_in("${src}" rev-parse HEAD)
set(first "${git_out}")
build_and_expect("first commit, built again" "${first}" false)

# What went wrong before: commits that change no source, and a build directory that is built again, not configured.
file(APPEND "${src}/notes.txt" "two\n")
git_in("${src}" commit -q -a -m second)
git_in("${src}" rev-parse HEAD)
set(second "${git_out}")
if("${second}" STREQUAL "${first}")
  give_up("the second commit has the hash of the first")
endif()
build_and_expect("second commit, built again" "${second}" false)
if(build_ok AND "${cost_not_checked}" STREQUAL "")
  string(REGEX MATCHALL "Building CXX object [^\n]*" compiled "${build_out}")
  set(more "")
  foreach(line IN LISTS compiled)
    if(NOT line MATCHES "git_stamp\\.cpp")
      string(APPEND more "${line}\n")
    endif()
  endforeach()
  if(NOT "${more}" STREQUAL "")
    check_failed("a new commit must cost one small compile (git_stamp.cpp) and a relink of each program, but the build also compiled" "${more}")
  else()
    message(STATUS "ok: only the stamp was compiled")
  endif()
endif()

build_and_expect("nothing changed, built again" "${second}" false)
if(build_ok AND "${cost_not_checked}" STREQUAL "")
  if(build_out MATCHES "Building CXX object|Linking CXX")
    check_failed("a build with nothing new must not compile or link anything, but it printed:" "${build_out}")
  else()
    message(STATUS "ok: nothing was compiled or linked")
  endif()
endif()
if(NOT "${cost_not_checked}" STREQUAL "")
  message(STATUS "not checked: what a rebuild compiles and links (${cost_not_checked})")
endif()

file(APPEND "${src}/notes.txt" "three\n")
build_and_expect("tracked file modified, built again" "${second}" true)

file(REMOVE_RECURSE "${work}")
