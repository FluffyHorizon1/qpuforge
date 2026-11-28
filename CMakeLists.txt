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
foreach(t qpuforge_cli:src/qpuforge.cpp test_qpuforge:tests/test_qpuforge.cpp design_study:examples/design_study.cpp alpha1_review:examples/alpha1_review.cpp alpha1_v11_review:examples/alpha1_v11_review.cpp beta1_review:examples/beta1_review.cpp gamma1_review:examples/gamma1_review.cpp decoder_study:examples/decoder_study.cpp decoder_sweeps:examples/decoder_sweeps.cpp surgery_study:examples/surgery_study.cpp test_fp_contract:tests/test_fp_contract.cpp)
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
# layout_presets reads the committed cell files (designs/*/layout_cells.json)
target_compile_definitions(test_qpuforge PRIVATE QPUFORGE_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
# `qpuforge --version` prints this (the version comes from project(... VERSION ...) above)
target_compile_definitions(qpuforge_cli PRIVATE QPUFORGE_VERSION="${PROJECT_VERSION}")
foreach(name alpha1_v11_review beta1_review gamma1_review decoder_study decoder_sweeps surgery_study)
  target_compile_definitions(${name} PRIVATE QPUFORGE_GIT_SHA="${QPUFORGE_GIT_SHA}" QPUFORGE_GIT_DIRTY=${QPUFORGE_GIT_DIRTY})
endforeach()
# the materials group reads its forge.materials.v0 sidecar fixtures from the source tree
target_compile_definitions(test_qpuforge PRIVATE QPUFORGE_TEST_DATA="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")

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
qpuforge_cli_test(veto_needs_bursts "qpuforge qec: --veto needs --bursts" qec ${CLI_CHIP} --veto)
qpuforge_cli_test(veto_options_need_veto "qpuforge qec: --impact-p needs --veto" qec ${CLI_CHIP} --bursts --impact-p 0.1)
qpuforge_cli_test(herald_options_need_herald "qpuforge qec: --herald-rule needs --herald" qec ${CLI_CHIP} --herald-rule detector)
qpuforge_cli_test(window_needs_stream "qpuforge qec: --window needs --decoder uf-stream" qec ${CLI_CHIP} --window 4)
qpuforge_cli_test(commit_below_window "qpuforge qec: --commit \\(4\\) must be smaller than --window \\(4\\)"
  qec ${CLI_CHIP} --decoder uf-stream --window 4 --commit 4)
qpuforge_cli_test(surgery_readout_needs_cnot "qpuforge surgery: --readout applies to --op cnot only" surgery --readout Z --shots 0)
qpuforge_cli_test(surgery_even_distance "qpuforge surgery: --d expects odd distances from 3 to 15, got 4" surgery --d 3,4 --shots 0)
qpuforge_cli_test(unknown_option_surgery "qpuforge surgery: unknown option '--shot' \\(did you mean '--shots'\\?\\)" surgery --shot 100)
qpuforge_cli_test(bad_number_surgery_shots "qpuforge surgery: --shots expects an integer >= 0, got '1e5'" surgery --shots 1e5)
qpuforge_cli_test(bad_seed_surgery "qpuforge surgery: --seed expects an integer from 0 to 18446744073709551615, got '-1'" surgery --shots 0 --seed -1)
qpuforge_cli_test(bad_choice_surgery "qpuforge surgery: --op expects zz, xx or cnot, got 'zzz'" surgery --op zzz)
qpuforge_cli_test(surgery_extra_argument "qpuforge surgery: unexpected argument 'extra'" surgery extra)
qpuforge_cli_test(bad_probability "qpuforge qec: --herald-fidelity expects a number from 0 to 1, got '1.5'" qec ${CLI_CHIP} --herald --herald-fidelity 1.5)
qpuforge_cli_test(loss_tangent_needs_materials "qpuforge analyze: --loss-tangent needs --materials <sidecar.json>"   # before the chip is read
  analyze ${CLI_NO_CHIP} --loss-tangent)
# accepted: one run per command with every option given, three for qec (the burst options; the decoding options;
# the streaming decoder);
# cz and gds stop at the missing chip file
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
# the Phase Gamma decoding options (union-find, heralds, impact veto), none at its default: the run echoes each value
qpuforge_cli_test(ok_qec_decoding
  "decoder: weighted union-find.*fidelity 0\\.900, false-herald rate 2\\.0e-03; detector rule, factor 0\\.50.*impact veto, memory-Z \\(2 rounds\\): window 4 rounds x radius 1\\.5, .*target 2\\.0e-02\\).*impacts injected in 10\\.0% of shots \\(onset uniform in the shot, local\\)"
  qec ${CLI_CHIP} --rounds 2 --shots 50 --seed 3 --decoder uf --herald --herald-fidelity 0.9 --herald-false 0.002
  --herald-factor 0.5 --herald-rule detector --bursts --veto --veto-window 4 --veto-radius 1.5 --veto-false 0.02
  --impact-p 0.1 --burst-radius 2)
qpuforge_cli_test(ok_qec_stream "decoder: streaming union-find.*window 4 rounds, committing 2 per step.*logical error per round"
  qec ${CLI_CHIP} --rounds 4 --shots 50 --seed 3 --decoder uf-stream --window 4 --commit 2)
qpuforge_cli_test(ok_surgery "CNOT d=3 .*noiseless tableau check, 2 random branches: 0 detector firings, 0 wrong observables, 0 wrong reference bits  \\(ok\\)"
  surgery --op cnot --d 3 --input 00 --readout Z --chip beta1 --shots 0 --seed 1 --merge-rounds 0 --pre 1 --post 1 --samples 2
  --cache ${CLI_OUT}_unused_cache.json)   # --shots 0: the noiseless check only (no CZ calibration)
# the ZZ and XX measurements (ok_surgery is the CNOT); a "FAILED" or an error in the output is a regression
qpuforge_cli_test(ok_surgery_zz "ZZ measurement d=3 \\|00>.*0 detector firings, 0 wrong observables, 0 wrong reference bits  \\(ok\\)"
  surgery --op zz --d 3 --input 00 --shots 0 --seed 5)
qpuforge_cli_test(ok_surgery_xx "XX measurement d=3 \\|\\+\\+>.*0 detector firings, 0 wrong observables, 0 wrong reference bits  \\(ok\\)"
  surgery --op xx --d 3 --shots 0)
set_tests_properties(cli_ok_surgery cli_ok_surgery_zz cli_ok_surgery_xx PROPERTIES FAIL_REGULAR_EXPRESSION "FAILED|error:")
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
# output files that cannot be written are errors (exit 2), not "wrote ..."
qpuforge_cli_test(write_fails_new "error: cannot write .*no-such-dir/chip.json" new ${CLI_OUT}_no-such-dir/chip.json --line 2)
qpuforge_cli_test(write_fails_layout "error: cannot write .*no-such-dir/chip.svg" layout ${CLI_OUT}_line2.json ${CLI_OUT}_no-such-dir/chip.svg)
qpuforge_cli_test(write_fails_hwgate_out "error: cannot write .*no-such-dir/gate.json"
  hwgate ${CLI_CHIP} ${CLI_OUT}_synth.json --shots 50 --out ${CLI_OUT}_no-such-dir/gate.json)
set_tests_properties(cli_write_fails_layout PROPERTIES FIXTURES_REQUIRED cli_line2)
set_tests_properties(cli_write_fails_hwgate_out PROPERTIES FIXTURES_REQUIRED cli_synth)
if(EXISTS /dev/full)   # opens, then every write fails: caught by the check after close()
  qpuforge_cli_test(write_fails_full "error: cannot write /dev/full" synth-chip ${CLI_CHIP} /dev/full --shots 50 --rounds 2 --z-only)
endif()
add_test(NAME cli_exit_status COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:qpuforge_cli> -DCHIP=${CLI_CHIP} -DMISSING=${CLI_NO_CHIP}
         -DOUT=${CLI_OUT}_exit_status_synth.json
         -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_cli_exit_status.cmake)
set_tests_properties(cli_exit_status PROPERTIES LABELS unit)
# usage on request (exit 0, stdout): the pass-regex tests cannot see the exit status, so a script checks it
foreach(m long:cli_help short:cli_help_short word:cli_help_word command:cli_help_command)
  string(REPLACE ":" ";" mp ${m})
  list(GET mp 0 mode)
  list(GET mp 1 tname)
  add_test(NAME ${tname} COMMAND ${CMAKE_COMMAND} -DEXE=$<TARGET_FILE:qpuforge_cli> -DMODE=${mode}
           -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/check_cli_help.cmake WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(${tname} PROPERTIES LABELS unit)
endforeach()
# CLI guards for `analyze --materials` (not validation checks), on the committed Alpha-1 v1.1 design: a sidecar with
# values is applied before the device report, atomforge's placeholder sidecar changes nothing, the loss tangent needs
# --loss-tangent, a chip design given as the sidecar is refused, and the usage text lists the flags.
set(QPUFORGE_CLI_CHIP ${CMAKE_CURRENT_SOURCE_DIR}/designs/alpha1-v11/alpha1_v11.json)
set(QPUFORGE_CLI_DATA ${CMAKE_CURRENT_SOURCE_DIR}/tests/data)
add_test(NAME cli_materials COMMAND qpuforge_cli analyze ${QPUFORGE_CLI_CHIP} --materials ${QPUFORGE_CLI_DATA}/materials-values.json)
add_test(NAME cli_materials_placeholder COMMAND qpuforge_cli analyze ${QPUFORGE_CLI_CHIP} --materials ${QPUFORGE_CLI_DATA}/materials-contract-CuAu.json)
add_test(NAME cli_materials_loss_tangent COMMAND qpuforge_cli analyze ${QPUFORGE_CLI_CHIP} --materials ${QPUFORGE_CLI_DATA}/materials-loss-tangent.json --loss-tangent)
add_test(NAME cli_materials_refused COMMAND qpuforge_cli analyze ${QPUFORGE_CLI_CHIP} --materials ${QPUFORGE_CLI_CHIP})
add_test(NAME cli_materials_loss_tangent_opt_in COMMAND qpuforge_cli analyze ${QPUFORGE_CLI_CHIP} --materials ${QPUFORGE_CLI_DATA}/materials-loss-tangent.json)
add_test(NAME cli_usage_materials COMMAND qpuforge_cli)
set_tests_properties(cli_materials PROPERTIES FAIL_REGULAR_EXPRESSION "error:" PASS_REGULAR_EXPRESSION
  "applied +Q_dielectric +materials.Q_dielectric 3e\\+06 -> 1e\\+06.*applied +junction_gap_ueV +materials.junction_gap_ueV 180 -> 190.*materials: Q_diel=1e\\+06  junction gap 190 ueV.*mean T1")
set_tests_properties(cli_materials_placeholder PROPERTIES FAIL_REGULAR_EXPRESSION "applied |error:" PASS_REGULAR_EXPRESSION
  "skipped +Q_dielectric +placeholder.*skipped +loss_tangent +placeholder.*skipped +junction_gap_ueV +placeholder.*materials: Q_diel=3e\\+06  junction gap 180 ueV.*mean T1")
set_tests_properties(cli_materials_loss_tangent PROPERTIES FAIL_REGULAR_EXPRESSION "error:" PASS_REGULAR_EXPRESSION
  "applied +loss_tangent +materials.Q_dielectric 3e\\+06 -> 4e\\+06, Q_dielectric = 1/loss_tangent.*materials: Q_diel=4e\\+06.*mean T1")
set_tests_properties(cli_materials_loss_tangent_opt_in PROPERTIES FAIL_REGULAR_EXPRESSION "applied +loss_tangent|error:" PASS_REGULAR_EXPRESSION
  "skipped +loss_tangent .*is opt-in.*materials: Q_diel=3e\\+06.*mean T1")
set_tests_properties(cli_materials_refused PROPERTIES PASS_REGULAR_EXPRESSION "error: .*alpha1_v11.json: materials sidecar: schema must be")
set_tests_properties(cli_usage_materials PROPERTIES PASS_REGULAR_EXPRESSION "analyze <chip.json> \\[--materials <sidecar.json> \\[--loss-tangent\\]\\]")
set_tests_properties(cli_materials cli_materials_placeholder cli_materials_loss_tangent cli_materials_loss_tangent_opt_in cli_materials_refused
  cli_usage_materials PROPERTIES LABELS unit)
# --materials on the other design-consuming sub-commands (it used to be silently ignored there): applied before the command
# runs, the sidecar's error is raised, and layout/gds refuse it instead of ignoring it.
add_test(NAME cli_materials_fab COMMAND qpuforge_cli fab ${QPUFORGE_CLI_CHIP} --materials ${QPUFORGE_CLI_DATA}/materials-values.json --samples 5)
add_test(NAME cli_materials_run COMMAND qpuforge_cli run ${QPUFORGE_CLI_CHIP} ghz:3 --traj 10 --shots 2 --materials ${QPUFORGE_CLI_DATA}/materials-values.json)
add_test(NAME cli_materials_budget COMMAND qpuforge_cli budget ${QPUFORGE_CLI_CHIP} ghz:3 --traj 10 --materials ${QPUFORGE_CLI_DATA}/materials-values.json)
add_test(NAME cli_materials_run_refused COMMAND qpuforge_cli run ${QPUFORGE_CLI_CHIP} ghz:3 --traj 10 --shots 2 --materials ${QPUFORGE_CLI_CHIP})
add_test(NAME cli_materials_layout_refused COMMAND qpuforge_cli layout ${QPUFORGE_CLI_CHIP} cli_materials_layout.svg --materials ${QPUFORGE_CLI_DATA}/materials-values.json)
add_test(NAME cli_loss_tangent_alone COMMAND qpuforge_cli fab ${QPUFORGE_CLI_CHIP} --loss-tangent --samples 5)
add_test(NAME cli_materials_gds_refused COMMAND qpuforge_cli gds ${QPUFORGE_CLI_CHIP} cli_materials_gds.gds --materials ${QPUFORGE_CLI_DATA}/materials-values.json)
add_test(NAME cli_materials_qec COMMAND qpuforge_cli qec ${CMAKE_CURRENT_SOURCE_DIR}/examples/alpha1.json --rounds 2 --shots 50 --seed 3 --materials ${QPUFORGE_CLI_DATA}/materials-values.json --loss-tangent)
add_test(NAME cli_materials_cz_missing_chip COMMAND qpuforge_cli cz ${CMAKE_CURRENT_BINARY_DIR}/cli_no_chip.json --materials ${QPUFORGE_CLI_DATA}/materials-values.json)
add_test(NAME cli_materials_hwgate_missing_chip COMMAND qpuforge_cli hwgate ${CMAKE_CURRENT_BINARY_DIR}/cli_no_chip.json ${QPUFORGE_CLI_DATA}/materials-values.json --materials ${QPUFORGE_CLI_DATA}/materials-values.json)
set_tests_properties(cli_materials_fab PROPERTIES FAIL_REGULAR_EXPRESSION "error:" PASS_REGULAR_EXPRESSION
  "applied +Q_dielectric +materials.Q_dielectric 3e\\+06 -> 1e\\+06.*fabrication Monte Carlo: 5 chips")
set_tests_properties(cli_materials_run PROPERTIES FAIL_REGULAR_EXPRESSION "error:" PASS_REGULAR_EXPRESSION
  "applied +Q_dielectric +materials.Q_dielectric 3e\\+06 -> 1e\\+06")
set_tests_properties(cli_materials_budget PROPERTIES FAIL_REGULAR_EXPRESSION "error:" PASS_REGULAR_EXPRESSION
  "applied +junction_gap_ueV.*scale with the gap.*error budget for GHZ-3")
set_tests_properties(cli_materials_run_refused PROPERTIES PASS_REGULAR_EXPRESSION "error: .*alpha1_v11.json: materials sidecar: schema must be")
set_tests_properties(cli_materials_layout_refused PROPERTIES PASS_REGULAR_EXPRESSION "qpuforge layout: --materials is not supported by `layout`")
set_tests_properties(cli_loss_tangent_alone PROPERTIES PASS_REGULAR_EXPRESSION "qpuforge fab: --loss-tangent needs --materials <sidecar.json>")
set_tests_properties(cli_materials_fab cli_materials_run cli_materials_budget cli_materials_run_refused cli_materials_layout_refused cli_loss_tangent_alone PROPERTIES LABELS unit)
set_tests_properties(cli_materials_gds_refused PROPERTIES PASS_REGULAR_EXPRESSION "qpuforge gds: --materials is not supported by `gds`")
# examples/alpha1.json has no physical coupler: no CZ calibration, so the run takes seconds
set_tests_properties(cli_materials_qec PROPERTIES FAIL_REGULAR_EXPRESSION "error:" PASS_REGULAR_EXPRESSION
  "applied +Q_dielectric +materials.Q_dielectric 3e\\+06 -> 1e\\+06.*memory-Z experiment, distance 3")
# cz and hwgate are expensive (CZ calibration, a full gate run): a chip file that does not exist shows that --materials was accepted by the
# command line and read the chip first (the sidecar is applied after the chip has loaded, so "cannot open" comes first)
set_tests_properties(cli_materials_cz_missing_chip cli_materials_hwgate_missing_chip PROPERTIES PASS_REGULAR_EXPRESSION "error: cannot open .*cli_no_chip.json")
set_tests_properties(cli_materials_gds_refused cli_materials_qec cli_materials_cz_missing_chip cli_materials_hwgate_missing_chip PROPERTIES LABELS unit)
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
  "leakage_frame:qec" "calibrated_noise_plumbing:qec" "layout:layout" "layout_geometry:layout" "layout_presets:layout" "json_rfc:unit" "json_depth:unit" "gds_write:layout"
  "matching:qec" "fault_enumeration:qec" "memory_bases:qec" "bursts:qec" "characterisation:qec" "designs_v04:physics"
  "union_find:qec" "heralds:qec" "impact_veto:qec" "streaming_uf:qec" "surgery:qec" "surgery_cnot:qec")
# fast unit groups for the serialisation layers (JSON parser/writer + design files, GDSII stream writer)
list(APPEND QPUFORGE_TESTS "json:unit" "gds:unit")
# forge.materials.v0 sidecar reader (atomforge -> qpuforge, docs/DATA_CONTRACT.md); fixtures in tests/data
list(APPEND QPUFORGE_TESTS "materials:unit")
# q18 opt-in flux-line distortion term (docs/design/flux-line-distortion.md); reads a copy of results/cz_cache.json
list(APPEND QPUFORGE_TESTS "flux_line:pulse")
# q17 FEM reference cross-check scaffold (fem.hpp): Elmer / Palace importers; synthetic fixtures in tests/data/fem
list(APPEND QPUFORGE_TESTS "fem:layout")
foreach(t ${QPUFORGE_TESTS})
  string(REPLACE ":" ";" parts ${t})
  list(GET parts 0 name)
  list(GET parts 1 label)
  add_test(NAME ${name} COMMAND test_qpuforge ${name})
  set_tests_properties(${name} PROPERTIES LABELS ${label} FAIL_REGULAR_EXPRESSION "FAIL ")
endforeach()
# End-to-end smoke runs of the review programs: 200 shots, outputs in the build tree, committed files untouched.
# They reuse results/cz_cache.json (copied); with a stale or missing cache they re-calibrate (a few minutes).
foreach(name alpha1_v11_review beta1_review gamma1_review surgery_study)
  add_test(NAME smoke_${name} COMMAND ${name} --smoke --out ${CMAKE_CURRENT_BINARY_DIR}/smoke_${name}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(smoke_${name} PROPERTIES LABELS smoke TIMEOUT 1800)
endforeach()
