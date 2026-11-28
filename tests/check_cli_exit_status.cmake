# CLI guard (not one of the validation checks): the exit statuses that the usage text and README document. The other
# cli_* tests match output with PASS_REGULAR_EXPRESSION, which ignores the exit status.
#   0 done, 1 usage error, 2 error while running, 3 hwgate: gate not passed (surgery: a noiseless check failed)
# Option dependencies and value ranges are checked on a chip file that does not exist: a usage error exits 1 before any
# file is read, so a check that went missing (or ran after load_design) shows as exit 2, "cannot open".
# Usage: cmake -DEXE=<qpuforge> -DCHIP=<examples/alpha1.json> -DMISSING=<a path that does not exist>
#              -DOUT=<a characterisation file to write> -P check_cli_exit_status.cmake
foreach(v EXE CHIP MISSING OUT)
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
# 3: a virtual chip with 8 us of dephasing the model does not know fails the gate by a wide margin (logical error per
# round above the limit, measured/model ratio about 4-5 against 2). Made here, so no committed design file is involved.
expect_status(0 synth-chip "${CHIP}" "${OUT}" --seed 1 --shots 200 --rounds 2,4,6 --z-only --hidden-dephasing-us 8)
expect_status(3 hwgate "${CHIP}" "${OUT}" --shots 200)

# qec: options that need another option
expect_status(1 qec "${MISSING}" --veto)
foreach(o --burst-rate --burst-tau-ms --burst-t1-us --burst-tphi-us)
  expect_status(1 qec "${MISSING}" ${o} 1)          # needs --bursts
endforeach()
foreach(o --veto-window --veto-radius --veto-false --impact-p --burst-radius)
  expect_status(1 qec "${MISSING}" --bursts ${o} 1)   # needs --veto
endforeach()
foreach(o --herald-fidelity --herald-false --herald-factor)
  expect_status(1 qec "${MISSING}" ${o} 0.5)        # needs --herald
endforeach()
expect_status(1 qec "${MISSING}" --herald-rule detector)
# qec: value ranges (probabilities and the herald factor in [0, 1], the veto window >= 1, radii and rates >= 0)
foreach(o --herald-fidelity --herald-false --herald-factor)
  expect_status(1 qec "${MISSING}" --herald ${o} 1.5)
  expect_status(1 qec "${MISSING}" --herald ${o} -0.1)
endforeach()
foreach(o --veto-false --impact-p)
  expect_status(1 qec "${MISSING}" --bursts --veto ${o} 1.5)
  expect_status(1 qec "${MISSING}" --bursts --veto ${o} -0.1)
endforeach()
expect_status(1 qec "${MISSING}" --bursts --veto --veto-window 0)
foreach(o --veto-radius --burst-radius)
  expect_status(1 qec "${MISSING}" --bursts --veto ${o} -1)
endforeach()
foreach(o --burst-rate --burst-tau-ms)
  expect_status(1 qec "${MISSING}" --bursts ${o} -1)
endforeach()
expect_status(1 qec "${MISSING}" --herald --herald-rule detecter)   # unknown choice (was exit 2)
# qec: the streaming decoder's window options need it, and it commits fewer rounds than it holds
expect_status(1 qec "${MISSING}" --window 4)
expect_status(1 qec "${MISSING}" --commit 2)
expect_status(1 qec "${MISSING}" --decoder uf-stream --window 3 --commit 3)
expect_status(1 qec "${MISSING}" --decoder uf-stream --window 1)
expect_status(1 qec "${MISSING}" --decoder uf-stream --commit 0)
# surgery (no chip file): the noiseless check alone passes; bad values are usage errors
expect_status(0 surgery --op cnot --d 3 --shots 0 --samples 1)
expect_status(1 surgery --op yy)
expect_status(1 surgery --d 3,4)
expect_status(1 surgery --readout X)                 # zz reads each patch in its input basis
expect_status(1 surgery --post 0)
expect_status(0 surgery --op zz --d 3 --gap 3 --decoder joint --shots 0 --samples 1)   # routing strip
expect_status(1 surgery --gap 2)                     # the strip length is odd
expect_status(1 surgery --gap 0)
expect_status(1 surgery --gap 17)
expect_status(1 surgery --decoder exact)             # qec's decoder names are not surgery's
# analyze: the materials options
expect_status(1 analyze "${MISSING}" --loss-tangent)
expect_status(1 analyze "${MISSING}" --materials "")
