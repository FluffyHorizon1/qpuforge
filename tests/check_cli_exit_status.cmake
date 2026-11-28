# CLI guard (not one of the validation checks): the exit statuses that the usage text and README document. The other
# cli_* tests match output with PASS_REGULAR_EXPRESSION, which ignores the exit status.
#   0 done, 1 usage error, 2 error while running (3, hwgate: gate not passed, needs a full gate run and is not checked)
# Usage: cmake -DEXE=<qpuforge> -DCHIP=<examples/alpha1.json> -DMISSING=<a path that does not exist> -P check_cli_exit_status.cmake
foreach(v EXE CHIP MISSING)
  if(NOT DEFINED ${v})
    message(FATAL_ERROR "check_cli_exit_status.cmake: -D${v}=... is required")
  endif()
endforeach()

function(expect_status status)
  execute_process(COMMAND "${EXE}" ${ARGN} RESULT_VARIABLE rc OUTPUT_VARIABLE out ERROR_VARIABLE err)
  if(NOT "${rc}" STREQUAL "${status}")
    message(SEND_ERROR "qpuforge ${ARGN}: exit status ${rc}, expected ${status}\n${err}")
  else()
    message(STATUS "exit ${rc}: qpuforge ${ARGN}")
  endif()
endfunction()

expect_status(0 --version)
expect_status(0 analyze "${CHIP}")
expect_status(1)                                    # no arguments: the usage
expect_status(1 anaylze "${CHIP}")                  # unknown command
expect_status(1 analyze "${CHIP}" --shots 5)        # unknown option
expect_status(1 run "${CHIP}" ghz:3 --shots)        # missing value
expect_status(1 run "${CHIP}" ghz:3 --shots lots)   # malformed number (was 2, from std::stoi)
expect_status(1 run "${CHIP}")                      # missing <circuit>
expect_status(2 analyze "${MISSING}")               # error while running: the chip file cannot be read
