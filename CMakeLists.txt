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
# CLI guards for `analyze --materials` (not validation checks), on the committed Alpha-1 v1.1 design: a sidecar with
# values is applied before the device report, atomforge's placeholder sidecar changes nothing, the loss tangent needs
# --loss-tangent, a chip design given as the sidecar is refused, and the usage text lists the flags.
set(QPUFORGE_CLI_CHIP ${CMAKE_CURRENT_SOURCE_DIR}/designs/alpha1-v11/alpha1_v11.json)
set(QPUFORGE_CLI_DATA ${CMAKE_CURRENT_SOURCE_DIR}/tests/data)
add_test(NAME cli_materials COMMAND qpuforge_cli analyze ${QPUFORGE_CLI_CHIP} --materials ${QPUFORGE_CLI_DATA}/materials-values.json)
add_test(NAME cli_materials_placeholder COMMAND qpuforge_cli analyze ${QPUFORGE_CLI_CHIP} --materials ${QPUFORGE_CLI_DATA}/materials-contract-CuAu.json)
add_test(NAME cli_materials_loss_tangent COMMAND qpuforge_cli analyze ${QPUFORGE_CLI_CHIP} --materials ${QPUFORGE_CLI_DATA}/materials-loss-tangent.json --loss-tangent)
add_test(NAME cli_materials_refused COMMAND qpuforge_cli analyze ${QPUFORGE_CLI_CHIP} --materials ${QPUFORGE_CLI_CHIP})
add_test(NAME cli_usage_materials COMMAND qpuforge_cli)
set_tests_properties(cli_materials PROPERTIES FAIL_REGULAR_EXPRESSION "error:" PASS_REGULAR_EXPRESSION
  "applied +Q_dielectric +materials.Q_dielectric 3e\\+06 -> 1e\\+06.*applied +junction_gap_ueV +materials.junction_gap_ueV 180 -> 190.*materials: Q_diel=1e\\+06  junction gap 190 ueV.*mean T1")
set_tests_properties(cli_materials_placeholder PROPERTIES FAIL_REGULAR_EXPRESSION "applied |error:" PASS_REGULAR_EXPRESSION
  "skipped +Q_dielectric +placeholder.*skipped +loss_tangent +placeholder.*skipped +junction_gap_ueV +placeholder.*materials: Q_diel=3e\\+06  junction gap 180 ueV.*mean T1")
set_tests_properties(cli_materials_loss_tangent PROPERTIES FAIL_REGULAR_EXPRESSION "error:" PASS_REGULAR_EXPRESSION
  "applied +loss_tangent +materials.Q_dielectric 3e\\+06 -> 4e\\+06, Q_dielectric = 1/loss_tangent.*materials: Q_diel=4e\\+06.*mean T1")
set_tests_properties(cli_materials_refused PROPERTIES PASS_REGULAR_EXPRESSION "error: .*alpha1_v11.json: materials sidecar: schema must be")
set_tests_properties(cli_usage_materials PROPERTIES PASS_REGULAR_EXPRESSION "analyze <chip.json> \\[--materials <sidecar.json> \\[--loss-tangent\\]\\]")
set_tests_properties(cli_materials cli_materials_placeholder cli_materials_loss_tangent cli_materials_refused cli_usage_materials PROPERTIES LABELS unit)
# CLI guards for `surgery` (not validation checks): --shots 0 runs the noiseless tableau check only (seconds, no Monte Carlo);
# malformed or misplaced values are errors (exit 2, ignored by a pass regex) and the usage text documents the exit statuses.
add_test(NAME cli_surgery_zz COMMAND qpuforge_cli surgery --op zz --d 3 --shots 0)
add_test(NAME cli_surgery_cnot COMMAND qpuforge_cli surgery --op cnot --d 3 --input 00 --shots 0)
add_test(NAME cli_surgery_bad_shots COMMAND qpuforge_cli surgery --shots 1e5)
add_test(NAME cli_surgery_bad_seed COMMAND qpuforge_cli surgery --shots 0 --seed -1)
add_test(NAME cli_surgery_bad_distance COMMAND qpuforge_cli surgery --d 4 --shots 0)
add_test(NAME cli_surgery_bad_readout COMMAND qpuforge_cli surgery --op zz --readout Z --shots 0)
add_test(NAME cli_usage_exit_status COMMAND qpuforge_cli)
set_tests_properties(cli_surgery_zz PROPERTIES FAIL_REGULAR_EXPRESSION "FAILED|error:" PASS_REGULAR_EXPRESSION
  "ZZ measurement d=3 \\|00>.*noiseless tableau check, 4 random branches: 0 detector firings, 0 wrong observables, 0 wrong reference bits  \\(ok\\)")
set_tests_properties(cli_surgery_cnot PROPERTIES FAIL_REGULAR_EXPRESSION "FAILED|error:" PASS_REGULAR_EXPRESSION
  "noiseless tableau check, 4 random branches: 0 detector firings, 0 wrong observables, 0 wrong reference bits  \\(ok\\)")
set_tests_properties(cli_surgery_bad_shots PROPERTIES PASS_REGULAR_EXPRESSION "error: --shots: '1e5' is not an integer")
set_tests_properties(cli_surgery_bad_seed PROPERTIES PASS_REGULAR_EXPRESSION "error: --seed: '-1' is not a non-negative integer")
set_tests_properties(cli_surgery_bad_distance PROPERTIES PASS_REGULAR_EXPRESSION "error: --d: ")
set_tests_properties(cli_surgery_bad_readout PROPERTIES PASS_REGULAR_EXPRESSION "error: --readout applies to --op cnot only")
set_tests_properties(cli_usage_exit_status PROPERTIES PASS_REGULAR_EXPRESSION "exit status: 0 done; 1 usage error;.*3 hwgate: gate not passed, surgery: the noiseless tableau check failed")
set_tests_properties(cli_surgery_zz cli_surgery_cnot cli_surgery_bad_shots cli_surgery_bad_seed cli_surgery_bad_distance cli_surgery_bad_readout cli_usage_exit_status PROPERTIES LABELS unit)
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
