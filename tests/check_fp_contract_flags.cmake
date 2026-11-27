# q13: the FMA guard must share compile flags with the real targets.
# If test_fp_contract lists -ffp-contract=off privately, it cannot detect the flag
# being dropped from qpuforge / the review programs / test_qpuforge.
if(NOT SOURCE_DIR)
  message(FATAL_ERROR "SOURCE_DIR is not set")
endif()
file(READ "${SOURCE_DIR}/CMakeLists.txt" txt)
if(NOT txt MATCHES [[set\(QPUFORGE_COMPILE_OPTIONS]])
  message(FATAL_ERROR "CMakeLists.txt must define QPUFORGE_COMPILE_OPTIONS shared by every target")
endif()
if(NOT txt MATCHES [[QPUFORGE_COMPILE_OPTIONS[^\n]*-ffp-contract=off]])
  message(FATAL_ERROR "QPUFORGE_COMPILE_OPTIONS must include -ffp-contract=off")
endif()
if(txt MATCHES [[target_compile_options\(test_fp_contract[^\n]*-ffp-contract=]])
  message(FATAL_ERROR "test_fp_contract has its own -ffp-contract flag; use QPUFORGE_COMPILE_OPTIONS so the guard tracks the real targets")
endif()
if(NOT txt MATCHES [[target_compile_options\(\$\{name\} PRIVATE \$\{QPUFORGE_COMPILE_OPTIONS\}\)]]
    AND NOT txt MATCHES [[test_fp_contract:tests/test_fp_contract.cpp]])
  message(FATAL_ERROR "the target helper must use QPUFORGE_COMPILE_OPTIONS (or test_fp_contract must be built by that helper)")
endif()
if(NOT txt MATCHES [[test_fp_contract:tests/test_fp_contract.cpp]]
    AND NOT txt MATCHES [[target_compile_options\(test_fp_contract PRIVATE \$\{QPUFORGE_COMPILE_OPTIONS\}\)]])
  message(FATAL_ERROR "test_fp_contract is not built with QPUFORGE_COMPILE_OPTIONS")
endif()
