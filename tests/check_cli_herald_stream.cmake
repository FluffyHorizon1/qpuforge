# CLI guard (not one of the validation checks): `qpuforge qec --decoder uf-stream --herald` prints the logical error
# "decoded without heralds" next to the run's own. Both must come from the streaming decoder. Before the fix the
# comparison copied the options that carry no streaming decoder and ran the batch one.
#   1. heralds that never fire (fidelity 0, false rate 0): the run and its comparison decode the same syndromes with the
#      same decoder, so the two P_L are equal;
#   2. that P_L differs from the batch union-find decoder's on the same seed (window 2 / commit 1 is far from batch), so
#      check 1 can tell the two decoders apart.
# Usage: cmake -DEXE=<qpuforge> -DCHIP=<examples/alpha1.json> -P check_cli_herald_stream.cmake
foreach(v EXE CHIP)
  if(NOT DEFINED ${v})
    message(FATAL_ERROR "check_cli_herald_stream.cmake: -D${v}=... is required")
  endif()
endforeach()
set(COMMON qec "${CHIP}" --rounds 6 --shots 4000 --seed 7 --herald --herald-fidelity 0 --herald-false 0)
execute_process(COMMAND "${EXE}" ${COMMON} --decoder uf-stream --window 2 --commit 1 RESULT_VARIABLE rc OUTPUT_VARIABLE stream ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "qpuforge qec --decoder uf-stream --herald: exit status ${rc}\n${err}")
endif()
execute_process(COMMAND "${EXE}" ${COMMON} --decoder uf RESULT_VARIABLE rc OUTPUT_VARIABLE batch ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "qpuforge qec --decoder uf --herald: exit status ${rc}\n${err}")
endif()
# the result row: "      6        860     0.0451    0.02725 +- 0.00257 ..." (rounds, round_ns, detection, P_L)
function(row_pl text out)
  if(NOT text MATCHES "\n +6 +[0-9]+ +[0-9.]+ +([0-9.]+) \\+- ")
    message(FATAL_ERROR "no result row for 6 rounds in:\n${text}")
  endif()
  set(${out} "${CMAKE_MATCH_1}" PARENT_SCOPE)
endfunction()
row_pl("${stream}" pl_stream)
row_pl("${batch}" pl_batch)
if(NOT stream MATCHES "decoded without heralds \\(streaming union-find\\): P_L ([0-9.]+),")
  message(FATAL_ERROR "the streaming run does not say that its comparison used the streaming decoder:\n${stream}")
endif()
set(pl_blind "${CMAKE_MATCH_1}")
message(STATUS "streaming: P_L ${pl_stream}, decoded without heralds ${pl_blind}; batch union-find: P_L ${pl_batch}")
if(NOT pl_blind STREQUAL pl_stream)
  message(FATAL_ERROR "heralds never fire, yet the comparison run gives P_L ${pl_blind} and the streaming run ${pl_stream}: they used different decoders")
endif()
if(pl_stream STREQUAL pl_batch)
  message(FATAL_ERROR "streaming and batch union-find agree (P_L ${pl_stream}): this guard cannot tell them apart; change its window or seed")
endif()
