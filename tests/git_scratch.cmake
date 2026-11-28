# Shared by the two source-stamp guards, tests/check_git_stamp.cmake and tests/check_git_stamp_wiring.cmake. Each sets
# `work`, its scratch directory in the build tree, and then includes this file, which
#   - empties ${work} and creates ${work}/home;
#   - keeps git away from the user's and the system's configuration (identity, commit signing, hooks, templates) and
#     from another repository named by the environment, and gives the commits a fixed author;
#   - sets `git` to the git the test uses: GIT_EXECUTABLE if that file exists, otherwise the git on PATH, otherwise
#     empty. This is the choice cmake/git_stamp.cmake makes;
#   - defines the functions below.
# The environment set here is inherited by everything the test starts: git, the stamp script, a build.
if("${work}" STREQUAL "")
  message(FATAL_ERROR "git_scratch.cmake: set work before including this file")
endif()
file(REMOVE_RECURSE "${work}")
file(MAKE_DIRECTORY "${work}/home")
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

set(git "${GIT_EXECUTABLE}")
if(NOT git OR NOT EXISTS "${git}")
  find_program(GIT_SCRATCH_GIT NAMES git PATHS ENV PATH NO_DEFAULT_PATH)
  set(git "${GIT_SCRATCH_GIT}")
endif()

# A failed check: <text>, then, if given, the output it is about, as it was printed. The script goes on to the
# remaining checks, and the test fails.
function(check_failed text)
  set_property(GLOBAL PROPERTY GIT_SCRATCH_CHECK_FAILED TRUE)
  if(ARGC GREATER 1)
    string(REPLACE "\n" "\n  " output "${ARGV1}")   # message() does not reflow an indented line
    string(APPEND text "\n  ${output}")
  endif()
  message(SEND_ERROR "${text}")
endfunction()

# A failure of the scaffolding (not of a check): removes the scratch directory, then stops. <why>, then the output.
function(give_up why)
  file(REMOVE_RECURSE "${work}")
  if(ARGC GREATER 1)
    string(REPLACE "\n" "\n  " output "${ARGV1}")
    string(APPEND why "\n  ${output}")
  endif()
  message(FATAL_ERROR "${why}")
endfunction()

# git <arguments> in <dir>; its output is left in git_out.
function(git_in dir)
  execute_process(COMMAND "${git}" ${ARGN} WORKING_DIRECTORY "${dir}"
                  RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err OUTPUT_STRIP_TRAILING_WHITESPACE)
  if(NOT "${rc}" STREQUAL "0")
    string(REPLACE ";" " " arguments "${ARGN}")
    give_up("git ${arguments} failed in ${dir} (${rc})" "${out}\n${err}")
  endif()
  set(git_out "${out}" PARENT_SCOPE)
endfunction()

# The end of a test that cannot go on without git: removes the scratch directory and prints <line>, the line that
# makes CTest report the test as skipped (SKIP_REGULAR_EXPRESSION in CMakeLists.txt). CTest lets that line win over
# the exit status, so it is not printed if a check has failed: a check that fails without git must fail the test.
function(skip_without_git line)
  file(REMOVE_RECURSE "${work}")
  get_property(failed GLOBAL PROPERTY GIT_SCRATCH_CHECK_FAILED)
  if(NOT failed)
    message(STATUS "${line}")
  endif()
endfunction()
