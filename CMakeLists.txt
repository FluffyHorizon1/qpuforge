cmake_minimum_required(VERSION 3.16)
project(qpuforge VERSION 0.3.0 LANGUAGES CXX)
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
foreach(t qpuforge_cli:src/qpuforge.cpp test_qpuforge:tests/test_qpuforge.cpp design_study:examples/design_study.cpp alpha1_review:examples/alpha1_review.cpp alpha1_v11_review:examples/alpha1_v11_review.cpp beta1_review:examples/beta1_review.cpp)
  string(REPLACE ":" ";" parts ${t})
  list(GET parts 0 name)
  list(GET parts 1 src)
  add_executable(${name} ${src})
  target_link_libraries(${name} PRIVATE qpuforge)
  # -ffp-contract=off: no fused multiply-add contraction, so results do not depend on -march / -mfma (or on an
  # FMA-by-default target such as aarch64). With contraction, -march=native moves the Alpha-1 Gen-0 noise model
  # enough to change eps from 0.01200 to 0.01223 at 8 threads (see results/PROVENANCE.md).
  target_compile_options(${name} PRIVATE -Wall -Wextra -O3 -ffp-contract=off)
  if(NOT OpenMP_CXX_FOUND)
    # serial build: the '#pragma omp' lines are intentionally ignored
    target_compile_options(${name} PRIVATE -Wno-unknown-pragmas)
  endif()
endforeach()
set_target_properties(qpuforge_cli PROPERTIES OUTPUT_NAME qpuforge)
# `qpuforge --version` prints this (the version comes from project(... VERSION ...) above)
target_compile_definitions(qpuforge_cli PRIVATE QPUFORGE_VERSION="${PROJECT_VERSION}")
foreach(name alpha1_v11_review beta1_review)
  target_compile_definitions(${name} PRIVATE QPUFORGE_GIT_SHA="${QPUFORGE_GIT_SHA}" QPUFORGE_GIT_DIRTY=${QPUFORGE_GIT_DIRTY})
endforeach()

# One CTest entry per validation test, labelled by layer:  ctest -L unit | physics | pulse | qec | layout | smoke
#   ctest -LE smoke    everything except the end-to-end review smoke runs
enable_testing()
# Build-flag guard (not one of the validation checks): fails if this build contracts a*b+c into an FMA.
add_executable(test_fp_contract tests/test_fp_contract.cpp)
target_compile_options(test_fp_contract PRIVATE -Wall -Wextra -O3 -ffp-contract=off)
add_test(NAME fp_contract COMMAND test_fp_contract)
set_tests_properties(fp_contract PROPERTIES LABELS unit)
set(QPUFORGE_TESTS
  "decomposition:unit" "transmon:physics" "coupler_and_flux:physics" "zz_and_inverse_design:physics" "qasm:unit" "compiler:unit" "decoherence:physics"
  "exchange:physics" "trajectories_vs_exact:physics" "ghz_metric:physics" "xeb:physics" "dynamic_circuits:unit"
  "surface_geometry:qec" "qec_circuits:qec" "frame_vs_statevector:qec" "eigensolver:unit" "pulse_cz:pulse"
  "leakage_frame:qec" "calibrated_noise_plumbing:qec" "layout:layout" "json_depth:unit")
foreach(t ${QPUFORGE_TESTS})
  string(REPLACE ":" ";" parts ${t})
  list(GET parts 0 name)
  list(GET parts 1 label)
  add_test(NAME ${name} COMMAND test_qpuforge ${name})
  set_tests_properties(${name} PROPERTIES LABELS ${label} FAIL_REGULAR_EXPRESSION "FAIL ")
endforeach()
# End-to-end smoke runs of the review programs: 200 shots, outputs in the build tree, committed files untouched.
# They reuse results/cz_cache.json (copied); with a stale or missing cache they re-calibrate (a few minutes).
foreach(name alpha1_v11_review beta1_review)
  add_test(NAME smoke_${name} COMMAND ${name} --smoke --out ${CMAKE_CURRENT_BINARY_DIR}/smoke_${name}
           WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  set_tests_properties(smoke_${name} PROPERTIES LABELS smoke TIMEOUT 1800)
endforeach()
