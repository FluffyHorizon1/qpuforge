# CLI guard (not one of the validation checks): usage on request. `--help`, `-h` and `help` print the usage to stdout and
# exit 0; `<command> --help` / `-h` (and `help <command>`) print that command's usage only, to stdout, exit 0, with
# nothing on stderr. A pass regex ignores the exit status, so this script checks it.
# Usage: cmake -DEXE=<qpuforge> -DMODE=<long|short|word|command> -P check_cli_help.cmake
foreach(v EXE MODE)
  if(NOT DEFINED ${v})
    message(FATAL_ERROR "check_cli_help.cmake: -D${v}=... is required")
  endif()
endforeach()

# expect_usage(<expected stdout regex> <argument>...): exit 0, empty stderr, stdout matches
function(expect_usage regex)
  execute_process(COMMAND "${EXE}" ${ARGN} RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT "${rc}" STREQUAL "0")
    message(SEND_ERROR "qpuforge ${ARGN}: exit status ${rc}, expected 0\n${err}")
  elseif(NOT "${err}" STREQUAL "")
    message(SEND_ERROR "qpuforge ${ARGN}: stderr is not empty:\n${err}")
  elseif(NOT "${out}" MATCHES "${regex}")
    message(SEND_ERROR "qpuforge ${ARGN}: stdout does not match '${regex}':\n${out}")
  else()
    message(STATUS "exit 0, usage on stdout: qpuforge ${ARGN}")
  endif()
endfunction()

# the full usage: the banner, every command and the exit-status note
set(FULL "^qpuforge — superconducting QPU chip design[^\n]*\n.*qpuforge new .*qpuforge analyze .*qpuforge surgery .*qpuforge --version \\| -V.*exit status: 0 done; 1 usage error")

if(MODE STREQUAL "long")
  expect_usage("${FULL}" --help)
elseif(MODE STREQUAL "short")
  expect_usage("${FULL}" -h)
elseif(MODE STREQUAL "word")
  expect_usage("${FULL}" help)
  expect_usage("^usage:\n  qpuforge analyze " help analyze)
elseif(MODE STREQUAL "command")
  # one command's usage: starts with "usage:", names the command, and no other command's synopsis
  expect_usage("^usage:\n  qpuforge analyze <chip.json> \\[--materials <sidecar.json> \\[--loss-tangent\\]\\]" analyze --help)
  expect_usage("^usage:\n  qpuforge analyze <chip.json> " analyze -h)
  expect_usage("^usage:\n  qpuforge analyze <chip.json> " analyze examples/alpha1.json --help)   # no chip file is read
  expect_usage("^usage:\n  qpuforge qec <chip.json> .*--decoder exact\\|uf\\|uf-stream" qec --help)
  expect_usage("^usage:\n  qpuforge surgery \\[--op zz\\|xx\\|cnot\\]" surgery -h)
  expect_usage("^usage:\n  qpuforge new <chip.json> " new --help)
  expect_usage("^usage:\n  qpuforge run <chip.json> <circuit> " run --help)
  expect_usage("^usage:\n  qpuforge hwgate <chip.json> <data.json> " hwgate -h)
  expect_usage("^usage:\n  qpuforge gds <chip.json> <out.gds> " gds --help)
  # a value that is not a help request: `--help` after an option that takes a value is that option's (bad) value
  execute_process(COMMAND "${EXE}" run --off --help RESULT_VARIABLE rc OUTPUT_QUIET ERROR_QUIET)
  if("${rc}" STREQUAL "0")
    message(SEND_ERROR "qpuforge run --off --help: exit status 0, expected a usage error")
  endif()
  # the usage of one command does not list the others
  execute_process(COMMAND "${EXE}" analyze --help OUTPUT_VARIABLE out)
  if("${out}" MATCHES "qpuforge (new|run|qec|surgery) ")
    message(SEND_ERROR "qpuforge analyze --help lists other commands:\n${out}")
  endif()
else()
  message(FATAL_ERROR "check_cli_help.cmake: unknown MODE '${MODE}'")
endif()
