cmake_minimum_required(VERSION 3.16)
project(qpuforge VERSION 0.4.0 LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
find_package(OpenMP)

# provenance for the review programs (examples/review_io.hpp): source commit at configure time
find_package(Git QUIET)
set(QPUFORGE_GIT_SHA "unknown")
set(QPUFORGE_GIT_DIRTY 0)
if(GIT_FOUND)
  execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  OUTPUT_VARIABLE QPUFORGE_GIT_SHA OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE _git_rc)
  if(NOT _git_rc EQUAL 0)
    set(QPUFORGE_GIT_SHA "unknown")
  endif()
  execute_process(COMMAND ${GIT_EXECUTABLE} status --porcelain --untracked-files=no WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                  OUTPUT_VARIABLE _git_status OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(NOT "${_git_status}" STREQUAL "")
    set(QPUFORGE_GIT_DIRTY 1)
  endif()
endif()
add_library(qpuforge INTERFACE)
target_include_directories(qpuforge INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
if(OpenMP_CXX_FOUND)
  target_link_libraries(qpuforge INTERFACE OpenMP::OpenMP_CXX)
endif()
# -ffp-contract=off: no fused multiply-add contraction, so results do not depend on -march / -mfma (or on an
# FMA-by-default target such as aarch64). With contraction, -march=native moves the Alpha-1 Gen-0 noise model
# enough to change eps from 0.01200 to 0.01223 at 8 threads (see results/PROVENANCE.md).
# Every target, including the fp_contract guard, must use this list (see tests/check_fp_contract_flags.cmake).
set(QPUFORGE_COMPILE_OPTIONS -Wall -Wextra -O3 -ffp-contract=off)
foreach(t qpuforge_cli:src/qpuforge.cpp test_qpuforge:tests/test_qpuforge.cpp design_study:examples/design_study.cpp alpha1_review:examples/alpha1_review.cpp alpha1_v11_review:examples/alpha1_v11_review.cpp beta1_review:examples/beta1_review.cpp gamma1_review:examples/gamma1_review.cpp test_fp_contract:tests/test_fp_contract.cpp)
  string(REPLACE ":" ";" parts ${t})
  list(GET parts 0 name)
  list(GET parts 1 src)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE qpuforge)
  # Shared with test_fp_contract so a lost -ffp-contract=off fails the guard as well as the real binaries.
  target_compile_options(${name} PRIVATE ${QPUFORGE_COMPILE_OPTIONS})
  if(NOT OpenMP_CXX_FOUND)
    # serial build: the '#pragma omp' lines are intentionally ignored
    target_compile_options(${name} PRIVATE -Wno-unknown-pragmas)
  endif()
endforeach()
set_target_properties(qpuforge_cli PROPERTIES OUTPUT_NAME qpuforge)
# `qpuforge --version` prints this (the version comes from project(... VERSION ...) above)
target_compile_definitions(qpuforge_cli PRIVATE QPUFORGE_VERSION="${PROJECT_VERSION}")
foreach(name alpha1_v11_review beta1_review gamma1_review)
  target_compile_definitions(${name} PRIVATE QPUFORGE_GIT_SHA="${QPUFORGE_GIT_SHA}" QPUFORGE_GIT_DIRTY=${QPUFORGE_GIT_DIRTY})
endforeach()

# One CTest entry per validation test, labelled by layer:  ctest -L unit | physics | pulse | qec | layout | smoke
#   ctest -LE smoke    everything except the end-to-end review smoke runs
enable_testing()
# Build-flag guard (not one of the validation checks): fails if this build contracts a*b+c into an FMA.
# Built with QPUFORGE_COMPILE_OPTIONS via the target helper above, so dropping the flag fails this test.
add_test(NAME fp_contract COMMAND test_fp_contract)
set_tests_properties(fp_contract PROPERTIES LABELS unit)
# CLI guards (not validation checks): `--version` / `-V` print the project version, and the usage text lists them.
add_test(NAME cli_version COMMAND qpuforge_cli --version)
add_test(NAME cli_version_short COMMAND qpuforge_cli -V)
set_tests_properties(cli_version cli_version_short PROPERTIES LABELS unit PASS_REGULAR_EXPRESSION "^qpuforge ${PROJECT_VERSION}")
add_test(NAME cli_usage COMMAND qpuforge_cli)   # no arguments: prints usage and exits 1 (ignored with a pass regex)
set_tests_properties(cli_usage PROPERTIES LABELS unit PASS_REGULAR_EXPRESSION "qpuforge --version \\| -V")
# CLI argument guards. Each command takes only the options it reads; anything else is a usage error that names the
# problem and exits 1 before any file is read. A pass regex ignores the exit status, so cli_exit_status checks the
# documented statuses. All fast: a 2-qubit line chip written into the build tree, examples/alpha1.json (no CZ
# calibration) or, for cz and gds (minutes of work), a chip file that does not exist, whose "error: cannot open" shows
# that the command line itself was accepted.
set(CLI_CHIP ${CMAKE_CURRENT_SOURCE_DIR}/examples/alpha1.json)
set(CLI_OUT ${CMAKE_CURRENT_BINARY_DIR}/cli)                # prefix of the files the CLI tests write
set(CLI_NO_CHIP ${CMAKE_CURRENT_BINARY_DIR}/cli_no_chip.json)   # never written
function(qpuforge_cli_test name regex)
  add_test(NAME cli_${name} COMMAND qpuforge_cli ${ARGN})
  set_tests_properties(cli_${name} PROPERTIES LABELS unit PASS_REGULAR_EXPRESSION "${regex}")
endfunction()
# rejected
qpuforge_cli_test(unknown_option_hwgate "qpuforge hwgate: unknown option '--hidden-dephasing-us' \\(an option of synth-chip\\)"
  hwgate ${CLI_CHIP} ${CMAKE_CURRENT_SOURCE_DIR}/designs/alpha1-v11/example_characterisation_synthetic.json --hidden-dephasing-us 8)
qpuforge_cli_test(unknown_option_run "qpuforge run: unknown option '--shot' \\(did you mean '--shots'\\?\\)" run ${CLI_CHIP} ghz:3 --shot 1000)
qpuforge_cli_test(unknown_option_qec "qpuforge qec: unknown option '--round' \\(did you mean '--rounds'\\?\\)" qec ${CLI_CHIP} --round 4)
qpuforge_cli_test(unknown_option_budget "qpuforge budget: unknown option '--shots'" budget ${CLI_CHIP} ghz:3 --shots 100)
qpuforge_cli_test(unknown_option_gds "qpuforge gds: unknown option '--synthesise' \\(did you mean '--synthesize'\\?\\)"
  gds ${CLI_CHIP} ${CLI_OUT}.gds --synthesise)
qpuforge_cli_test(unknown_short_option "qpuforge new: unknown option '-grid' \\(did you mean '--grid'\\?\\)" new ${CLI_OUT}.json -grid 3x3)
qpuforge_cli_test(unknown_command "qpuforge: unknown command 'anaylze' \\(did you mean 'analyze'\\?\\)" anaylze ${CLI_CHIP})
qpuforge_cli_test(missing_value "qpuforge run: option '--shots' requires a value" run ${CLI_CHIP} ghz:3 --shots)
qpuforge_cli_test(missing_value_before_option "qpuforge cz: option '--cache' requires a value" cz ${CLI_CHIP} --cache --out ${CLI_OUT}_cal.json)
qpuforge_cli_test(bad_number_shots "qpuforge run: --shots expects an integer >= 0, got 'lots'" run ${CLI_CHIP} ghz:3 --shots lots)
qpuforge_cli_test(bad_number_exponent "qpuforge qec: --shots expects an integer >= 1, got '1e5'" qec ${CLI_CHIP} --shots 1e5)
qpuforge_cli_test(bad_number_real "qpuforge new: --f-hi expects a number, got '5.4GHz'" new ${CLI_OUT}.json --line 2 --f-hi 5.4GHz)
qpuforge_cli_test(bad_count "qpuforge fab: --samples expects an integer >= 1, got '0'" fab ${CLI_CHIP} --samples 0)
qpuforge_cli_test(bad_choice "qpuforge new: --arch expects fixed or tunable, got 'tunabel' \\(did you mean 'tunable'\\?\\)"
  new ${CLI_OUT}.json --line 2 --arch tunabel)
qpuforge_cli_test(bad_noise_source "qpuforge run: --off: 'contrl' is not one of .*\\(did you mean 'control'\\?\\)" run ${CLI_CHIP} ghz:3 --off zz,contrl)
qpuforge_cli_test(bad_circuit "qpuforge run: circuit 'ghz:three': ghz:N needs a whole number N >= 1" run ${CLI_CHIP} ghz:three)
qpuforge_cli_test(repeated_option "qpuforge run: option '--traj' given more than once" run ${CLI_CHIP} ghz:3 --traj 2 --traj 3)
qpuforge_cli_test(extra_argument "qpuforge analyze: unexpected argument 'extra.json'" analyze ${CLI_CHIP} extra.json)
qpuforge_cli_test(missing_argument "qpuforge run: missing <circuit>" run ${CLI_CHIP})
qpuforge_cli_test(preset_takes_only_distance "qpuforge new: --name cannot be combined with --preset" new ${CLI_OUT}.json --preset alpha1 --name X)
qpuforge_cli_test(burst_options_need_bursts "qpuforge qec: --burst-rate needs --bursts" qec ${CLI_CHIP} --burst-rate 1)
# accepted: one run per command, every option given (cz and gds stop at the missing chip file)
qpuforge_cli_test(ok_new "wrote .*cli_line2.json \\(2 qubits, 1 couplers\\)"
  new ${CLI_OUT}_line2.json --line 2 --f-hi 5.4 --f-lo 4.8 --alpha -0.22 --g-MHz 4 --arch tunable --Q 3e6 --plan 8
  --stagger-MHz 60 --name cli-line2 --feedlines 1 --kappa-MHz 1 --chi-MHz 0 --purcell 30 --readout-ns 500 --reset-ns 200
  --reset-err 0.003 --err2q 0.002 --residual-zz-kHz 5 --active-reset)
qpuforge_cli_test(ok_new_preset "wrote .*cli_gen0.json \\(17 qubits, 24 couplers\\)" new ${CLI_OUT}_gen0.json --preset gen0 --distance 3)
qpuforge_cli_test(ok_analyze "cli-line2 .+ 2 qubits, 1 couplers" analyze ${CLI_OUT}_line2.json)
qpuforge_cli_test(ok_layout "wrote .*cli_line2.svg" layout ${CLI_OUT}_line2.json ${CLI_OUT}_line2.svg)
qpuforge_cli_test(ok_fab "collision-free yield: [0-9]+ / 3 " fab ${CLI_OUT}_line2.json --samples 3 --seed 2)
qpuforge_cli_test(ok_run "state fidelity .*\\(2 trajectories\\)" run ${CLI_OUT}_line2.json ghz:2 --traj 2 --shots 4 --seed 5 --off zz,control --top 2)
qpuforge_cli_test(ok_budget "error budget for GHZ-2" budget ${CLI_OUT}_line2.json ghz:2 --traj 2 --seed 5)
set_tests_properties(cli_ok_new PROPERTIES FIXTURES_SETUP cli_line2)
set_tests_properties(cli_ok_analyze cli_ok_layout cli_ok_fab cli_ok_run cli_ok_budget PROPERTIES FIXTURES_REQUIRED cli_line2)
qpuforge_cli_test(ok_qec "logical-error floor from impacts"
  qec ${CLI_CHIP} --rounds 2 --shots 50 --seed 3 --basis Z --statevector 0 --off leakage --decoder exact --leakage-removal none
  --cache ${CLI_OUT}_unused_cache.json --bursts --burst-rate 0.1 --burst-tau-ms 25 --burst-t1-us 1 --burst-tphi-us 0)
qpuforge_cli_test(ok_synth_chip "wrote .*cli_synth.json: 17 qubits"
  synth-chip ${CLI_CHIP} ${CLI_OUT}_synth.json --seed 1 --shots 50 --rounds 2 --z-only --hidden-dephasing-us 0 --cache ${CLI_OUT}_unused_cache.json)
qpuforge_cli_test(ok_hwgate "wrote .*cli_gate.json"   # exit 3 (gate not passed) is ignored by the pass regex
  hwgate ${CLI_CHIP} ${CLI_OUT}_synth.json --phase alpha --shots 50 --seed 1 --out ${CLI_OUT}_gate.json --cache ${CLI_OUT}_unused_cache.json)
set_tests_properties(cli_ok_synth_chip PROPERTIES FIXTURES_SETUP cli_synth)
set_tests_properties(cli_ok_hwgate PROPERTIES FIXTURES_REQUIRED cli_synth)
qpuforge_cli_test(ok_cz "error: cannot open .*cli_no_chip.json"
  cz ${CLI_NO_CHIP} --out ${CLI_OUT}_cal.json --cache ${CLI_OUT}_unused_cache.json --g-qc 160 --f-on 6.1 --rise 7)
qpuforge_cli_test(ok_gds "error: cannot open .*cli_no_chip.json"
  gds ${CLI_NO_CHIP} ${CLI_OUT}.gds --cells ${CLI_OUT}_cells.json --synthesize --from ${CLI_OUT}_cells0.json --cell TOP)
add_test(NAME cli_exit_status COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:qpuforge_cli> -DCHIP=${CLI_CHIP} -DMISSING=${CLI_NO_CHIP}
         -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_cli_exit_status.cmake)
set_tests_properties(cli_exit_status PROPERTIES LABELS unit)
add_test(NAME fp_contract_flags COMMAND ${CMAKE_COMMAND}
  -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}
  -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_fp_contract_flags.cmake)
set_tests_properties(fp_contract_flags PROPERTIES LABELS unit)
# Thread-count determinism guard (not one of the validation checks): Monte Carlo results must be bit-identical for any
# OpenMP thread count (docs/design/deterministic-mc.md). 'determinism' compares 1/2/3/4/8 threads in one process;
# 'determinism_env' runs separate processes under OMP_NUM_THREADS=1,2,4,8 and compares their digests.
add_executable(test_determinism tests/test_determinism.cpp)
target_link_libraries(test_determinism PRIVATE qpuforge)
target_compile_options(test_determinism PRIVATE ${QPUFORGE_COMPILE_OPTIONS})
if(NOT OpenMP_CXX_FOUND)
  target_compile_options(test_determinism PRIVATE -Wno-unknown-pragmas)
endif()
add_test(NAME determinism COMMAND test_determinism)
add_test(NAME determinism_env COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:test_determinism>
         -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/compare_thread_digests.cmake)
set_tests_properties(determinism determinism_env PROPERTIES LABELS unit)
set(QPUFORGE_TESTS
  "decomposition:unit" "transmon:physics" "coupler_and_flux:physics" "zz_and_inverse_design:physics" "qasm:unit" "compiler:unit" "decoherence:physics"
  "exchange:physics" "trajectories_vs_exact:physics" "ghz_metric:physics" "xeb:physics" "dynamic_circuits:unit"
  "surface_geometry:qec" "qec_circuits:qec" "frame_vs_statevector:qec" "eigensolver:unit" "pulse_cz:pulse"
  "leakage_frame:qec" "calibrated_noise_plumbing:qec" "layout:layout" "json_rfc:unit" "json_depth:unit" "gds_write:layout"
  "matching:qec" "fault_enumeration:qec" "memory_bases:qec" "bursts:qec" "characterisation:qec" "designs_v04:physics")
# fast unit groups for the serialisation layers (JSON parser/writer + design files, GDSII stream writer)
list(APPEND QPUFORGE_TESTS "json:unit" "gds:unit")
foreach(t ${QPUFORGE_TESTS})
  string(REPLACE ":" ";" parts ${t})
  list(GET parts 0 name)
  list(GET parts 1 label)
  add_test(NAME ${name} COMMAND test_qpuforge ${name})
  set_tests_properties(${name} PROPERTIES LABELS ${label} FAIL_REGULAR_EXPRESSION "FAIL ")
endforeach()
# End-to-end smoke runs of the review programs: 200 shots, outputs in the build tree, committed files untouched.
# They reuse results/cz_cache.json (copied); with a stale or missing cache they re-calibrate (a few minutes).
foreach(name alpha1_v11_review beta1_review gamma1_review)
  add_test(NAME smoke_${name} COMMAND ${name} --smoke --out ${CMAKE_CURRENT_BINARY_DIR}/smoke_${name}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(smoke_${name} PROPERTIES LABELS smoke TIMEOUT 1800)
endforeach()
